import CircuitElementGroup
import Bus


class BusGroup(CircuitElementGroup):
    dss_module_name, ele_name = 'Bus', 'Bus'

    def __init__(self, dss):
        super().__init__(self, dss)

    def _collect_elements(self, dss):
        super()._collect_elements(self, dss)
        # get phases for each bus from Node Names
        for node_name in dss.Circuit.AllNodeNames():
            name, phase = node_name.split('.')
            bus = self._name_to_object_dict[name]
            bus.phases.append(phase)

