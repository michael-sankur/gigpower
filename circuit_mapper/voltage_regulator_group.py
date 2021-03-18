from line_group import LineGroup
from voltage_regulator import VoltageRegulator
from typing import List


class VoltageRegulatorGroup(LineGroup):
    dss_module_name = 'RegControls'
    ele_class = VoltageRegulator

    def __init__(self, dss, line_group):
        super().__init__(dss, line_group=line_group)

    def _add_edge(self, vr):
        """
        Adds the upstream and downstream lines of voltage regulators to
        self.adj and self.reverse_adj
        """
        super()._add_edge(vr.upstream_line)
        super()._add_edge(vr.downstream_line)

    def get_tx_buses(self) -> List:
        """
        overwrite super to get tx buses of upstream, downstream lines
        in vr index order
        resulting list length is 2 * self.num_elements
        even number indices are upstream_line tx bus indices
        odd number indices are downstream_line tx bus indices
        """
        tx_buses = []
        for vr in self.get_elements():
            tx_buses.append(vr.upstream_line.tx)
            tx_buses.append(vr.downstream_line.tx)
        return tx_buses

    def get_rx_buses(self) -> List:
        """
        overwrite super to get rx buses of upstream, downstream lines
        in vr index order
        resulting list length is 2 * self.num_elements
        even number indices are upstream_line rx bus indices
        odd number indices are downstream_line rx bus indices
        """
        tx_buses = []
        for vr in self.get_elements():
            tx_buses.append(vr.upstream_line.rx)
            tx_buses.append(vr.downstream_line.rx)
        return tx_buses
