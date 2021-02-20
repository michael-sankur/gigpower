import CircuitElement
import numpy as np

class Line(CircuitElement):

    def __init__(self, name: str, txnode_name, rxnode_name):
        super().__init__(name)
        self.key = (txnode_name, rxnode_name)  # tuple of (txnode_name, rxnode_name)
        self.length = 0.0
        self.FZpu = np.zeros((3, 3), dtype='complex')
        self.Vbase = 0.0
        self.Ibase = 0.0
        self.Zbase = 0.0
        self.voltageRegulators = []  # hold a list of voltage regulators

    def __str__(self) -> str:
        return super().__str__() + f' txnode, rxnode: {self.key}'
