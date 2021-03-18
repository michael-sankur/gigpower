from circuit_element_group import CircuitElementGroup
from typing import Tuple, List
from line import Line
import numpy as np


class LineGroup(CircuitElementGroup):
    dss_module_name = 'Lines'
    ele_class = Line

    def __init__(self, dss, bus_group, **kwargs):
        """
        init self.adj, self.reverse_adj, self.buses
        then call super
        """
        self.adj = {}  # adjacency matrix -> { bus_name: [downstream buses]}
        # reverse adjacency matrix -> { bus_name: [upstream buses]}
        self.reverse_adj = {}
        self.buses = bus_group
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

    def add_element(self, line):
        """ Call super(), and add the new line's topology"""
        super().add_element(line)
        self._add_edge(line)

    def _add_edge(self, line):
        """ Adds line's topology to self.adj and self.reverse_adj"""
        tx_bus, rx_bus = line.key
        if tx_bus not in self.adj:
            self.adj[tx_bus] = [rx_bus]
        else:
            self.adj[tx_bus].append(rx_bus)
        if rx_bus not in self.reverse_adj:
            self.reverse_adj[rx_bus] = [tx_bus]
        else:
            self.reverse_adj[rx_bus].append(tx_bus)

    def get_bus_ph_matrix(self) -> np.ndarray:
        """
        5 x self.num_elements matrix, columns indexed by line index
        row 0: tx_bus_index
        row 1: rx_bus_index
        row 2 - 4: A, B, C phases
        """
        bp_matrix = np.zeros((self.num_elements, 5), dtype=int)
        for line in self.get_elements():
            bp_idx = self.get_idx(line.__name__)
            tx_bus_name, rx_bus_name = line.key
            bp_matrix[bp_idx][0:2] = np.asarray([self.buses.get_idx(
                tx_bus_name), self.buses.get_idx(rx_bus_name)])
            bp_matrix[bp_idx][2:] = line.get_phase_matrix()
        return bp_matrix.transpose()

    def get_num_lines_x_ph(self):
        """ returns sum of active phases across all lines """
        return sum([len(line.phases) for line in self.get_elements()])
