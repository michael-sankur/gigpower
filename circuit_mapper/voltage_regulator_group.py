from circuit_element_group import CircuitElementGroup

class VoltageRegulatorGroup(CircuitElementGroup):
    dss_module_name, ele_name = 'RegControls', 'voltage_regulator.VoltageRegulator'

    def __init__(self, dss, name, line_group):
        """
        Pass the LineGroup to the constructor so that
        Voltage Regulators can be assigned to Lines.
        """
        super().__init__(dss, name)
        self._collect_elements(dss, line_group)

    def _collect_elements(self, dss, line_group):
        """
        Assign VoltageRegulators to synthetic Lines in the line_group.
        Set Line phases based on Voltage Regulators.
        """
        super()._collect_elements(self, dss)  # collects VoltageRegulators
        for vreg in self.get_elements():
            syn_line = line_group.get_line_from_key(vreg.key)
            if not syn_line.voltageRegulators:
                syn_line.voltageRegulators = []
            # syn line phases are the union of all vreg phases
            syn_line.phases = list(set(syn_line.phases + vreg.phases))
            # add voltage regulator to this line
            syn_line.voltageRegulators.append(vreg)
