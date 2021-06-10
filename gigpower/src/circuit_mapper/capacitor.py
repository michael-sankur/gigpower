from . circuit_element import CircuitElement


class Capacitor(CircuitElement):
    dss_module_name = 'Capacitors'

    def __init__(self, name: str, dss):
        super().__init__(name, dss)
        cappu = dss.Capacitors.kvar() * 1000 / self.Sbase / len(self.phases)
        self._set_attr_val_by_phase('cappu', cappu)
