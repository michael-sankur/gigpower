import CircuitElementGroup


class CapacitorGroup(CircuitElementGroup):
    dss_module_name, ele_name = 'Capacitors', 'Capacitor'

    def __init__(self, dss):
        super().__init__(self, dss)

    def get_sum_cappu_by_bus(self):
        pass
