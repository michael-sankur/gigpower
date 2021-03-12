from circuit_element_group import CircuitElementGroup
from capacitor import Capacitor


class CapacitorGroup(CircuitElementGroup):
    dss_module_name = 'Capacitors'
    ele_class = Capacitor

    def get_sum_cappu_by_bus(self):
        pass
