from circuit_element import CircuitElement
import numpy as np


class Capacitor(CircuitElement):
    dss_module_name = 'Capacitors'

    def __init__(self, name: str, dss):
        super().__init__(self, name, dss)
        cappu = dss.Capcitor.kvar() * 1000 / self.Sbase / len(self.phases)
        self.cappu = np.zeros(3)
        self.cappu[np.asarray(self.phases) - 1] = cappu

