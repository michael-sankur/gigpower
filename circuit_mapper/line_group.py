from circuit_element_group import CircuitElementGroup
from voltage_regulator_group import VoltageRegulatorGroup
import numpy as np
from typing import Tuple
from line import Line


class LineGroup(CircuitElementGroup):
    dss_module_name = 'Lines'
    ele_class = Line

    def __init__(self, dss, *args):
        """
        Call CircuitElementGroup.__init__, then map adjacency matrix
        and map voltage regulators
        """
        super().__init__(dss, *args)
        self._key_to_element_dict = {line.key: line for line in self.get_elements()}
        self._set_topology()
        #  if this is a regular line group, and not any other subclass,
        # map the voltage regulators.
        if self.__class__ == 'LineGroup':
            self.voltage_regulators = VoltageRegulatorGroup(dss, self)

    def _set_topology(self):
        """ Computes adjacency matrix from (tx, rx) edges in the Line Group"""
        self.adj = {}  # adjacency matrix -> { bus_name: [downstream buses]}
        # reverse adjacency matrix -> { bus_name: [upstream buses]}
        self.reverse_adj = {}
        for line in self.get_elements():
            tx_bus, rx_bus = line.key
            if tx_bus not in self.adj:
                self.adj[tx_bus] = [rx_bus]
            else:
                self.adj[tx_bus].append(rx_bus)
            if rx_bus not in self.reverse_adj:
                self.reverse_adj[rx_bus] = [tx_bus]
            else:
                self.reverse_adj[rx_bus].append(tx_bus)

    def get_line_from_key(self, key: Tuple[str, str]):
        """ return the Line with the key (tx_bus, rx_bus)"""
        return self._key_to_element_dict[key]
