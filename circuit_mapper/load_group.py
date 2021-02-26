from circuit_element_group import CircuitElementGroup


class LoadGroup(CircuitElementGroup):
    dss_module_name, ele_name = 'Loads', 'load.Load'

    def __init__(self, dss):
        super().__init__(dss)

    def _collect_elements(self, dss):
        self._names = dss.Loads.AllNames()
        for name in self._names:
            load = Load(name, dss)
            self._name_to_object_dict[name] = load

    def get_sum_spu_by_bus(self):
        pass
