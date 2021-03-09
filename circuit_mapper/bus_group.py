from circuit_element_group import CircuitElementGroup
from bus import Bus


class BusGroup(CircuitElementGroup):
    dss_module_name = 'Bus'
    ele_class = Bus

    def __init__(self, dss):
        super().__init__(dss)

    def _collect_names(self, dss, **args):
        self._names = dss.Circuit.AllBusNames()
