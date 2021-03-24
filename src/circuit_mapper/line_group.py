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
        """ return the Line with the key (tx_bus, rx_bus)"""
        return self._key_to_element_dict[key]

    def get_bus_ids(self, which: str) -> List:
        """
        Return a list of tx or rx buses in this LineGroup, in line index order
        param which = 'tx' or 'rx'
        """
        return [getattr(line, which) for line in self.get_elements()]

    def add_element(self, line, **kwargs):
        """
        Call super(), and add the new line's topology.
        add line to self._key_to_element_dict
        note that adding a SyntheticLine will NOT increment self.num_elements
        """
        super().add_element(line, **kwargs)
        self._key_to_element_dict[line.key] = line
        self._add_edge(line)

    def _add_edge(self, line):
        """ Adds line's topology to self.adj and self.reverse_adj"""
        tx_bus, rx_bus = line.key
        try:
            if rx_bus not in self.adj[tx_bus]:
                self.adj[tx_bus].append(rx_bus)
        except KeyError:
            self.adj[tx_bus] = [rx_bus]

        try: 
            if tx_bus not in self.reverse_adj[rx_bus]:
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

    def get_idx(self, obj: Union[str, CircuitElement, Tuple]) -> int:
        """
        override super to check for tuples first
        """
        if isinstance(obj, tuple):
            return super().get_idx(self.get_line_from_key(obj))
        
        try: 
            return super().get_idx(obj)
        except KeyError:
            return super().get_idx(self.get_line_from_key(obj.key))

    def get_num_lines_x_ph(self):
        """ returns sum of active phases across all lines """
        return sum([len(line.phases) for line in self.get_elements()])

    def get_X_matrix(self):
        """ n x 9 matrix, indexed by Line index"""
        return self._get_attr_by_idx('xmat')

    def get_R_matrix(self):
        """ n x 9 matrix, indexed by Line index"""
        return self._get_attr_by_idx('rmat')

    def get_parents(self, bus_name: str) -> List[str]:
        """ returns a list of buses upstream of bus_name """
        if bus_name in self.reverse_adj:
            return self.reverse_adj[bus_name]
        return []

    def get_children(self, bus_name: str) -> List[str]:
        """ returns a list of buses downstream of bus_name """
        if bus_name in self.adj:
            return self.adj[bus_name]
        return []

  