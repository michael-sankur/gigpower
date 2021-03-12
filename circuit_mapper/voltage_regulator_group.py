from line_group import LineGroup
from line import SyntheticLine
from voltage_regulator import VoltageRegulator


class VoltageRegulatorGroup(LineGroup):
    dss_module_name = 'RegControls'
    ele_class = VoltageRegulator

    def __init__(self, dss, line_group):
        super().__init__(dss, line_group=line_group)

    def _collect_elements(self, dss, line_group):
        """
        Each voltage regulator is modeled with two lines:
        1. upstream: txBus --> regControlBus
        2. downstream: regControlBus--> rxBus
        opendss already has a Line for the downstream line, so it should
        be present in Circuit.lines
        Creates a SyntheticLine for the upstream line, find the downstream line
        from the line_group,  and assign voltage regulators to both lines
        """
        super()._collect_elements(dss)  # create VR, map self.adj and self.reverse_adj
        for vreg in self.get_elements():
            upstream_line = SyntheticLine(vreg)
            rx_line = self._find_downstream_line(vreg, line_group)
            downstream_line = SyntheticLine(rx_line)
            upstream_line.add_voltage_regulator(vreg)
            downstream_line.add_voltage_regulator(vreg)
            self._add_edge(upstream_line)
            self._add_edge(downstream_line)

    def _find_downstream_line(self, vreg, line_group):
        for line in line_group.get_elements():
            if line.tx == vreg.regControl_bus:
                return line
