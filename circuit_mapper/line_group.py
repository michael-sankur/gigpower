from circuit_element_group import CircuitElementGroup
from typing import Tuple, List
from line import Line


class LineGroup(CircuitElementGroup):
    dss_module_name = 'Lines'
    ele_class = Line

    def __init__(self, dss, **kwargs):
        """
        Call CircuitElementGroup.__init__, then map adjacency matrix
        and map voltage regulators
        """
        self.adj = {}  # adjacency matrix -> { bus_name: [downstream buses]}
        # reverse adjacency matrix -> { bus_name: [upstream buses]}
        self.reverse_adj = {}
        super().__init__(dss, **kwargs)

    def get_line_from_key(self, key: Tuple[str, str]):
        """ return the Line with the key (tx_bus, rx_bus)"""
        return self._key_to_element_dict[key]

    def get_tx_buses(self) -> List:
        """ Return a list of tx buses in this LineGroup, in line index order"""
        return [line.tx for line in self.get_elements()]

    def get_rx_buses(self) -> List:
        """ Return a list of rx buses in this LineGroup, in line index order"""
        return [line.rx for line in self.get_elements()]

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
