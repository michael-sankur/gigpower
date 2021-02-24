import CircuitElement
import numpy as np
from utils import pad_phases, parse_phases


class Line(CircuitElement):
    dss_module_name = 'Lines'

    def __init__(self, name: str, dss):
        super().__init__(name)
        self.length = dss.Lines.length()
        self._set_FZpu()

    def _set_related_bus(self, dss):
        """ override super class to save two buses"""
        dss.Lines.Name(self.__name__)
        self.tx, self.rx = dss.Lines.Bus1(), dss.Lines.Bus2()
        super()._set_base_vals(self.tx)  # set base vals based on tx bus
        super()._set_phases_from_bus(self.tx)  # set phases based on tx bus
        self.key = (self.tx, self.rx)

    def _set_FZpu(self, dss):
        fz_mult = 1 / self.Zbase * self.length  # use tx node's Z base
        num_phases = len(self.phases)
        RM = np.asarray(dss_data['RMatrix'])
        XM = np.asarray(dss_data['XMatrix'])
        ZM = fz_mult * (RM + 1j*XM)
        ZM = np.reshape(ZM, (ZM.shape[0]//num_phases, num_phases))  # reshape
        # pad the Z matrix
        self.FZpu = pad_phases(ZM, (3, 3), parse_phases(self.phases))

    def __str__(self) -> str:
        return super().__str__() + f' txnode, rxnode: {self.key}'
