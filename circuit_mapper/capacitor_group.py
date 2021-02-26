from circuit_element_group import CircuitElementGroup
from capacitor import Capacitor


class CapacitorGroup(CircuitElementGroup):
    dss_module_name = 'Capacitors'
    ele_class = Capacitor

    def __init__(self, dss):
        super().__init__(dss)

    def get_sum_cappu_by_bus(self):
        pass
