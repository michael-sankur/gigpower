import CircuitElementGroup
import VoltageRegulatorGroup
import numpy as np

class LineGroup(CircuitElementGroup):
    dss_module_name, ele_name = 'Lines', 'Line'

    def __init__(self, dss):
        """
        Call CircuitElementGroup.__init__, then map adjacency matrix
        and map voltage regulators
        """
        super().__init__(self, dss)
        self._key_to_element_dict = {line.key: line for line in self.get_elements()}
        self._set_topology()
        #  if this is a regular line group, and not any other subclass,
        # map the voltage regulators.
        if self.__class__ == 'LineGroup':
            self.voltage_regulators = VoltageRegulatorGroup(dss, self)

    def _set_topology(self):
        """ set adjacency matrices"""
        self.adj = {}  # adjacency matrix -> { bus_name: [downstream buses]}
        # reverse adjacency matrix -> { bus_name: [upstream buses]}
        self.reverse_adj = {}
        for line in self._name_to_idx_dict.values():
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

    def get_tx_idx_matrix(self, bus_group):
        """
        n x 1 matrix of tx bus indices. Indexed by line index,
        which is the same value as in opendss
        """
        idx_matrix = np.zeros(self.num_elements)
        for idx, line_name in self._idx_to_name_dict().items():
            tx_bus = self.get_element(line_name).tx
            idx_matrix[idx] = bus_group.get_idx(tx_bus)
        return idx_matrix

