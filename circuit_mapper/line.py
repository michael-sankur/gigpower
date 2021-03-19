from circuit_element import CircuitElement
import numpy as np
from utils import pad_phases, parse_phases, parse_dss_bus_name, \
    parse_dss_phases, parse_phase_matrix


class Line(CircuitElement):
    dss_module_name = 'Lines'

    def __init__(self, name: str, dss):
        super().__init__(name, dss)
        self._set_length(dss)
        self._set_FZpu(dss)
        self._set_x_r_matrices(dss)

    def _set_length(self, dss):
        dss.Lines.Name(self.__name__)
        self.length = dss.Lines.Length()

    def _set_related_bus(self, dss):
        """ override super class to save tx, rx, and two buses"""
        dss.Lines.Name(self.__name__)
        self.tx, self.rx = parse_dss_bus_name(dss.Lines.Bus1()), \
            parse_dss_bus_name(dss.Lines.Bus2())
        self.related_bus = self.tx
        self.key = (self.tx, self.rx)

    def _set_phases(self, dss):
        """ override super class to save phases from dss.Lines.Bus1()"""
        dss.Lines.Name(self.__name__)
        if dss.Lines.IsSwitch():  # assume 3-phase
            self.phases = ['1', '2', '3']
            self.phase_matrix = np.asarray([1, 1, 1])
        else:
            self.phases = parse_dss_phases(dss.Lines.Bus1())
            self.phase_matrix = parse_phase_matrix(self.phases)

    def _set_FZpu(self, dss):
        fz_mult = 1 / self.Zbase * self.length  # use tx node's Z base
        num_phases = len(self.phases)
        RM = np.asarray(dss.Lines.RMatrix())
        XM = np.asarray(dss.Lines.XMatrix())
        ZM = fz_mult * (RM + 1j*XM)
        ZM = np.reshape(ZM, (ZM.shape[0]//num_phases, num_phases))  # reshape
        # pad the Z matrix
        self.FZpu = pad_phases(ZM, (3, 3), parse_phases(self.phases))

    def add_voltage_regulator(self, vreg):
        try:
            self.voltage_regulators.append(vreg)
        except AttributeError:
            self.voltage_regulators = [vreg]

    def _set_x_r_matrices(self, dss):
        """ retrieve impedance and reactance matrices of a line
            pad by phase to 3x3 matrix, and flatten
            set self.xmat and self.rmat to 9x1 matrix
        """
        try:
            dss.Lines.Name(self.__name__)
            xmat = dss.Lines.XMatrix()
            rmat = dss.Lines.RMatrix()
        except Exception:  # for transformers, voltage regs, and Synthetic Lines
            xmat, rmat = np.zeros(9), np.zeros(9)
        line_phases = self.get_phase_matrix()
        if len(xmat) == 1:  # set the diagonals to the value
            self.xmat = np.identity(9) * xmat
            self.rmat = np.identity(9) * rmat
        elif len(xmat) == 9:
            self.xmat, self.rmat = xmat, rmat
        elif len(xmat) == 4:
            xmat = np.reshape(xmat, (2, 2))
            rmat = np.reshape(rmat, (2, 2))
            self.xmat = pad_phases(xmat, (3, 3), line_phases).flatten()
            self.rmat = pad_phases(rmat, (3, 3), line_phases).flatten()
        else:
            raise IndexError(f"Xmat is length {len(xmat)}, expected 1, 4, or 9 elements")

    def __str__(self) -> str:
        return super().__str__() + f' txnode, rxnode: {self.key}'


class SyntheticLine(Line):
    """
    Special class to handle voltage regulators and transformers
    Synthetic Lines are used to handle topology information about voltage regulators
    and transformers.
    """
    def __init__(self, ele):
        """ ele: a Transformer or Voltage Regulator"""
        self.__name__ = ele.__name__
        self.key = ele.key
        self.tx, self.rx = self.key
        self.related_bus = ele.related_bus
        self.length = 0
        self.phases = ele.phases
