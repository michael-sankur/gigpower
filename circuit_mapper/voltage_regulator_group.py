from line_group import LineGroup
from voltage_regulator import VoltageRegulator
from typing import List


class VoltageRegulatorGroup(LineGroup):
    dss_module_name = 'RegControls'
    ele_class = VoltageRegulator

    def __init__(self, dss, line_group):
        super().__init__(dss, bus_group=line_group.buses, line_group=line_group)

    def _add_edge(self, vr):
        """
        Adds the upstream and downstream lines of voltage regulators to
        self.adj and self.reverse_adj
        """
        super()._add_edge(vr.upstream_line)
        super()._add_edge(vr.downstream_line)

    def get_bus_ids(self, which: str) -> List:
        """
        overwrite super to get tx buses of upstream, downstream lines
        in vr index order
        param which = 'tx' or 'rx'
        resulting list length is 2 * self.num_elements
        even number indices are upstream_line tx bus indices
        odd number indices are downstream_line tx bus indices
        """
        buses = []
        for vr in self.get_elements():
            buses.append(getattr(vr.upstream_line, which))
            buses.append(getattr(vr.downstream_line, which))
        return buses
