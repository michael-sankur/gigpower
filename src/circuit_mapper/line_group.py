from . circuit_element_group import CircuitElementGroup
from . circuit_element import CircuitElement
from typing import Tuple, List, Union
from . line import Line
import numpy as np
from . utils import pad_phases


class LineGroup(CircuitElementGroup):
    dss_module_name = 'Lines'
    ele_class = Line

    def __init__(self, dss, bus_group, **kwargs):
        """
        init self.adj, self.reverse_adj, self.buses
        then call super
        """
        self.buses = bus_group
        # adjacency matrix -> { bus_name: [downstream buses]}
        self.adj = {}
        # reverse adjacency matrix -> { bus_name: [upstream buses]}
        self.reverse_adj = {}
        self._key_to_element_dict = {}
        super().__init__(dss, **kwargs)
  
    def get_line_from_key(self, key: Tuple[str, str]):
        """
        return the Line with the key (tx_bus, rx_bus)
        first searches self, then searches transformers, finally 
        voltage regulators
        """
        try:
            return self._key_to_element_dict[key]
        except KeyError:
            try:
                return self.transformers._key_to_element_dict[key]
            except KeyError:
                return self.voltage_regulators._key_to_element_dict[key]

    def get_bus_ids(self, which: str) -> List:
        """
        Return a list of tx or rx buses in this LineGroup, in line index order
        param which = 'tx' or 'rx'
        """
        return [getattr(line, which) for line in self.get_elements()]

    def add_element(self, line, unique_key=True):
        """
        Call super(), and add the new line's topology.
        add line to self._key_to_element_dict
        """
        super().add_element(line)
        if unique_key:
            self._key_to_element_dict[line.key] = line
        else:
            if line.key in self._key_to_element_dict:
                self._key_to_element_dict[line.key].append(line)
            else:
                self._key_to_element_dict[line.key] = [line]
        self._add_edge(line)

    def _add_edge(self, line):
        """ Adds line's topology to self.adj and self.reverse_adj"""
        tx_bus, rx_bus = line.key
        try:
            self.adj[tx_bus].append(rx_bus)
        except KeyError:
            self.adj[tx_bus] = [rx_bus]
        try:
            self.reverse_adj[rx_bus].append(tx_bus)
        except KeyError:
            self.reverse_adj[rx_bus] = [tx_bus]

    def get_bus_ph_matrix(self) -> np.ndarray:
        """
        5 x self.num_elements matrix, columns indexed by line index
        row 0: tx_bus_index
        row 1: rx_bus_index
        row 2 - 4: A, B, C phases
        """
        bp_matrix = np.zeros((self.num_elements, 5), dtype=int)
        for line in self.get_elements():
            bp_idx = self.get_idx(line)
            tx_bus_name, rx_bus_name = line.key
            bp_matrix[bp_idx][0:2] = np.asarray([self.buses.get_idx(
                tx_bus_name), self.buses.get_idx(rx_bus_name)])
            bp_matrix[bp_idx][2:] = line.phase_matrix
        return bp_matrix.transpose()

    def get_idx(self, obj: Union[str, CircuitElement, tuple]) -> int:
        """
        override super to include transformers and voltage regulators
        if this LineGroup has transformers and voltage regulators, will
        include them in indexing in the following order:
        [0, self.num_elements]: Lines
        [self.num_elements, transformers.num_elements -1 ]: Transformers
        [transformers.num_elements, voltage_regulators.num_elements-1]:
                                                VoltageRegulators
        """
        if isinstance(obj, tuple):
            if obj in self._key_to_element_dict:
                line = self._key_to_element_dict[obj]
                return self._name_to_idx_dict[line.__name__]
            # check transformers
            elif obj in self.transformers._key_to_element_dict:  
                return self.num_elements + self.transformers.get_idx(obj)
            # check vrs, will return a list of vrs
            elif obj in self.voltage_regulators._key_to_element_dict:  
                return self.voltage_regulators._key_to_element[obj]
        else:
            try:
                return super().get_idx(obj)  # search self
            except KeyError:  # check transforemrs
                try:
                    return self.num_elements + self.transformers.get_idx(obj)
                except KeyError:  # check vrs
                    return self.num_elements + self.transformers.num_elements \
                        + self.voltage_regulators.get_idx(obj)


    def get_num_lines_x_ph(self):
        """ returns sum of active phases across all lines """
        return sum([len(line.phases) for line in self.get_elements()])

    def get_X_matrix(self):
        """ n x 9 matrix, indexed by Line index"""
        return self._get_attr_by_idx('xmat', 'rows')

    def get_R_matrix(self):
        """ n x 9 matrix, indexed by Line index"""
        return self._get_attr_by_idx('rmat', 'rows')

    def get_upstream_buses(self, bus_name: str, inc_xfm=False,
                           inc_regs=False) -> List[str]:
        """ returns a list of buses upstream of bus_name """
        return self._get_adj(bus_name, 'upstream', inc_xfm=inc_xfm, inc_regs=inc_regs)

    def get_downstream_buses(self, bus_name: str, inc_xfm=False,
                             inc_regs=False) -> List[str]:
        """ returns a list of buses downstream of bus_name """
        return self._get_adj(bus_name, 'downstream', inc_xfm=inc_xfm, inc_regs = inc_regs)

    def _get_adj(self, bus_name: str, which, inc_xfm, inc_regs) -> List[str]:
        """ for get_children() and get_parents"""
        if which == 'upstream':
            target = 'reverse_adj'
        elif which == 'downstream':
            target = 'adj'

        self_adj = getattr(self, target)
        if bus_name in self_adj:
            return self_adj[bus_name]

        if inc_xfm:
            try:
                xfm_adj = getattr(self.transformers, target)
                if bus_name in self.transformers.reverse_adj:
                    return xfm_adj[bus_name]
            except AttributeError:
                pass

        if inc_regs:
            try:
                vr_adj = getattr(self.voltage_regulators, target)
                if bus_name in vr_adj:
                    return vr_adj[bus_name]
            except AttributeError:
                pass

        return []

    def get_line_list(self, bus_name: str, which, inc_xfm=False,
                      inc_regs=False) -> np.asarray:
        """
        list of line indices for lines going out or coming into bus_name
        param which: 'out' , returns lines starting with bus_name, 
        or 'in', returns lines terminating at bus_name
        """
        if which == 'out':
            buses = self.get_downstream_buses(bus_name, inc_xfm, inc_regs)
            if not buses:
                return []
            return [self.get_line_from_key((bus_name, rx)).__name__ for rx in buses]
        if which == 'in':
            buses = self.get_upstream_buses(bus_name, inc_xfm, inc_regs)
            if not buses:
                return []
            return [self.get_line_from_key((tx, bus_name)).__name__ for tx in buses]
