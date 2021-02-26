from circuit_element_group import CircuitElementGroup
from voltage_regulator import VoltageRegulator


class VoltageRegulatorGroup(CircuitElementGroup):
    dss_module_name = 'RegControls'
    ele_class = VoltageRegulator

    def __init__(self, dss, line_group):
        """
        Pass the LineGroup to the constructor so that
        Voltage Regulators can be assigned to Lines.
        """
        super().__init__(dss, line_group)

    def _collect_elements(self, dss, *args):
        """
        Assign VoltageRegulators to synthetic Lines in the line_group.
        Set Line phases based on Voltage Regulators.
        """
        super()._collect_elements(dss, *args)  # collects VoltageRegulators
        line_group = args[0]
        for vreg in self.get_elements():
            syn_line = line_group.get_line_from_key(vreg.key)
            if not syn_line.voltageRegulators:
                syn_line.voltageRegulators = []
            # syn line phases are the union of all vreg phases
            syn_line.phases = list(set(syn_line.phases + vreg.phases))
            # add voltage regulator to this line
            syn_line.voltageRegulators.append(vreg)
