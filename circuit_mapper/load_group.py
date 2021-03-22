from circuit_element_group import CircuitElementGroup
from load import Load
import numpy as np


class LoadGroup(CircuitElementGroup):
    dss_module_name, ele_class = 'Loads', Load
    # ZIP_V = np.asarray([0.10, 0.05, 0.85, 0.10, 0.05, 0.85, 0.80])
    # ZIP_V = np.asarray([0, 0, 1, 0, 0, 1, 0.75])  # constant power
    ZIP_V = np.asarray([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.75])  # constance impedance

    def __init__(self, dss, bus_group):
        self.buses = bus_group
        # set aPQ, aI, aZ

        # if dss.Loads.ZipV(): # TODO: uncomment after finished testing
        #     self.set_zip_values(dss.Loads.ZipV())

        # else:  # default to value on LoadGroup class
        self.set_zip_values(self.__class__.ZIP_V)
        super().__init__(dss)

    def load_powers(self, solution=None, zip_V: np.ndarray = None) -> np.ndarray:
        """
            returns load powers summed by Bus
            if solution and zip_V are not given, will default to nominal values

            param solution: optional Solution argument
            param zip_V: optional 6x1 zip values array

            returns: ndarray of shape (num_buses, 3), indexed by bus index
        """
        if not solution:
            bus_V_matrix = np.ones((self.buses.num_elements, 3), dtype=complex)
        else:
            pass
        load_powers = np.zeros((self.buses.num_elements, 3))

        for load in self.get_elements():
            bus_idx = self.buses.get_idx(load.related_bus)
            bus_V = bus_V_matrix[bus_idx]
            spu_real, spu_imag = load.ppu, load.qpu
            if not zip_V:
                aPQ_p, aI_p, aZ_p = load.aPQ_p, load.aI_p, load.aZ_p
                aPQ_q, aI_q, aZ_q = load.aPQ_q, load.aI_q, load.aZ_q
            else:
                aZ_p = zip_V[0]
                aI_p = zip_V[1]
                aPQ_p = zip_V[2]
                aZ_q = zip_V[3]
                aI_q = zip_V[4]
                aPQ_q = zip_V[5]

            for idx, ph in enumerate(load.phase_matrix):
                if ph == 1:
                    temp1 = aPQ_p + aI_p * abs(bus_V)[idx] + aZ_p * ((abs(bus_V))**2)[idx]
                    real = temp1 * spu_real

                    temp2 = aPQ_q + aI_q * abs(bus_V)[idx] + aZ_q * ((abs(bus_V))**2)[idx]
                    imag = temp2 * spu_imag

                    load_powers[idx] += real + (1j * imag)
        return load_powers

    def get_ppu_matrix(self):
        """
        return 3 x n matrix of load.ppu values summed over Bus
        columns indexed by Bus index, padded by phase
        """
        return self._get_attr_by_bus('ppu', orient='col')

    def get_qpu_matrix(self):
        """
        return 3 x n matrix of kvar values summed over Bus
        columns indexed by Bus index, padded by phase
        """
        return self._get_attr_by_bus('qpu', orient='col')

    def set_zip_values(self, zipV: np.ndarray):
        self.zipV = zipV
        # array mapping: [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min votlage pu]
        self.aZ_p, self.aI_p, self.aPQ_p = self.zipV[0:3]
        self.aZ_q, self.aI_q, self.aPQ_q = self.zipV[3:6]
