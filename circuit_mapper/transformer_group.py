from line_group import LineGroup


class TransformerGroup(LineGroup):
    dss_module_name, ele_name = 'Transformers', 'transformer.Transformer'

    def __init__(self, dss, line_group):
        """
        Pass the LineGroup to the constructor so that
        Transformers can be added to LineGroup
        """
        super().__init__(dss, line_group)

    def _collect_names(self, dss, line_group):
        """ Override LineGroup._collect_names to exclude voltage regulators"""
        transformers = dss.Transformers.AllNames()
        vregs = line_group.voltage_regulators.get_elements()
        # exclude the transformers that are paired with voltage regulators
        vr_transformers = [r.__name__ for r in vregs]
        self._names = [n for n in transformers if n not in vr_transformers]
