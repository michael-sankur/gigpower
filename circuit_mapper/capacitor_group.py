from circuit_element_group import CircuitElementGroup


class CapacitorGroup(CircuitElementGroup):
    dss_module_name, ele_name = 'Capacitors', 'capacitor.Capacitor'

    def __init__(self, dss):
        super().__init__(dss)

    def get_sum_cappu_by_bus(self):
        pass
