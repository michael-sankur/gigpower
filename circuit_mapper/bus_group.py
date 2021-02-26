from circuit_element_group import CircuitElementGroup
from bus import Bus


class BusGroup(CircuitElementGroup):
    dss_module_name = 'Bus'
    ele_class = Bus

    def __init__(self, dss):
        super().__init__(dss)

    def _collect_elements(self, dss, **args):
        for name in self._names:
            self._name_to_object_dict[name] = Bus(name, dss)

        # get phases for each bus from Node Names
        for node_name in dss.Circuit.AllNodeNames():
            name, phase = node_name.split('.')
            bus = self._name_to_object_dict[name]
            bus.phases += str(phase)

    def _collect_names(self, dss, **args):
        self._names = dss.Circuit.AllBusNames()
