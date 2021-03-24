from . circuit_element import CircuitElement
from typing import List
from . utils import parse_dss_bus_name, parse_dss_phases, parse_phase_matrix


class Load(CircuitElement):
    dss_module_name = 'Loads'

    def __init__(self, name, dss):
        super().__init__(name, dss)
        self.kW = dss.Loads.kW()
        self.kvar = dss.Loads.kvar()

        # divide ppu and qpu by number of phases
        # store as zero-padded 1x3 matrix
        ppu, qpu = self.kW / 1000 / len(self.phases), self.kvar / 1000 / len(self.phases)
        self._set_attr_val_by_phase('ppu', ppu)
        self._set_attr_val_by_phase('qpu', qpu)
        # set aPQ, aI, aZ
        if dss.Loads.ZipV():
            self.set_zip_values(dss.Loads.ZipV())
        else:  # default to constant impedance
            self.set_zip_values([1, 0, 0, 1, 0, 0])
        self.spu = self.ppu + 1j*self.qpu

    def __str__(self) -> str:
        return self.name

    def _set_related_bus(self, dss):
        dss.Loads.Name(self.__name__)
        self.related_bus = parse_dss_bus_name(dss.CktElement.BusNames()[0])

    def _set_phases(self, dss):
        dss.Loads.Name(self.__name__)
        self.phases = parse_dss_phases(dss.CktElement.BusNames()[0])
        self.phase_matrix = parse_phase_matrix(self.phases)

    def set_zip_values(self, zipV: List):
        self.zipV = zipV
        # array mapping: [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min votlage pu]
        self.aZ_p, self.aI_p, self.aPQ_p = self.zipV[0:3]
        self.aZ_q, self.aI_q, self.aPQ_q = self.zipV[3:6]

    def _set_kvar(self, kvar: float) -> None:
        """
        Reset a load's kvar and recalculate all load parameters based on new kvar.
        Note that for the bus associated with this load,
        Bus.sum_spu will not be updated by this method.
        Use Circuit.set_kvar() to also update the Bus.sum_spu.
        """
        self.kvar = kvar

        # divide by number of phases
        qpu = self.kvar / 1000 / len(self.phases)

        # set to 3x1 based on phases
        self._set_attr_val_by_phase('qpu', qpu)
        self.spu = self.ppu + 1j*self.qpu

    def _set_kW(self, kW: float) -> None:
        """
        Reset a load's kW and recalculate all load parameters based on new kvar.
        Note that for the bus associated with this load,
        Bus.sum_spu will not be updated by this method.
        Use Circuit.set_kvar() to also update the Bus.sum_spu.
        """
        self.kW = kW

        # divide by number of phases
        ppu = self.kW / 1000 / len(self.phases)

        self._set_attr_val_by_phase('ppu', ppu)
        self.spu = self.ppu + 1j*self.qpu
