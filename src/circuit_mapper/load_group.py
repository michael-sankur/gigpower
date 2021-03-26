from . circuit_element_group import CircuitElementGroup
from . load import Load
import numpy as np


class LoadGroup(CircuitElementGroup):
    dss_module_name, ele_class = 'Loads', Load

    def __init__(self, dss, bus_group):
        self.buses = bus_group
        super().__init__(dss)

    def get_powers(self, solution=None) -> np.ndarray:
        """
            returns load powers summed by Bus, given an optional voltage solution
            if solution is not given, will default to nominal values

            param solution: optional Solution argument

            returns: ndarray of shape (num_buses, 3), indexed by bus index
        """
        if not solution:
            bus_V_matrix = np.ones((self.buses.num_elements, 3), dtype=complex)
        else:
            bus_V_matrix = solution.V
        load_powers = np.zeros((self.buses.num_elements, 3))

        for load in self.get_elements():
            bus_idx = self.buses.get_idx(load.related_bus)
            bus_V = bus_V_matrix[bus_idx]
            spu_real, spu_imag = load.ppu, load.qpu

            aPQ_p, aI_p, aZ_p = self.aPQ_p, self.aI_p, self.aZ_p
            aPQ_q, aI_q, aZ_q = self.aPQ_q, self.aI_q, self.aZ_q

            for idx, ph in enumerate(load.phase_matrix):
                if ph == 1:
                    temp1 = aPQ_p + aI_p * abs(bus_V)[idx] + aZ_p * ((abs(bus_V))**2)[idx]
                    real = temp1 * spu_real

                    temp2 = aPQ_q + aI_q * abs(bus_V)[idx] + aZ_q * ((abs(bus_V))**2)[idx]
                    imag = temp2 * spu_imag

                    load_powers[idx] += real + (1j * imag)
        return load_powers

    def get_nominal_powers(self):
        return self.get_powers()

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

    # def _set_zip_values(self, ZIP_V):
    #     """
    #     Private helper method called by Circuit to set zip values
    #     as attributes on all LoadGroups
    #     according to what is set on the Circuit class.
    #     For the public, user-facing way to change zip values, use
    #     Solution.set_zip_values()
    #     """
    #     self.zipV = ZIP_V
    #     self.aZ_p, self.aI_p, self.aPQ_p = ZIP_V[0:3]
    #     self.aZ_q, self.aI_q, self.aPQ_q = ZIP_V[3:6]
    #     self.min_voltage_pu = ZIP_V[6]
