from circuit_element import CircuitElement
import numpy as np


class Capacitor(CircuitElement):
    dss_module_name = 'Capacitors'

    def __init__(self, name: str, dss):
        super().__init__(name, dss)
        cappu = dss.Capacitors.kvar() * 1000 / self.Sbase / len(self.phases)
        self.set_cappu(cappu)

    def set_cappu(self, cappu: float):
        """self.cappu to passed value, for phases in self.phases"""
        self.cappu = np.zeros(3)
        self.cappu[self.get_ph_idx_matrix()] = cappu
