from . circuit_element_group import CircuitElementGroup
from . load import Load
import numpy as np
import pandas as pd


class LoadGroup(CircuitElementGroup):
    dss_module_name, ele_class = 'Loads', Load

    def __init__(self, dss, bus_group, zip_v):
        self.buses = bus_group
        super().__init__(dss)
        self.ZIP_V = zip_v
        self.aZ_p, self.aI_p, self.aPQ_p = self.ZIP_V[0:3]
        self.aZ_q, self.aI_q, self.aPQ_q = self.ZIP_V[3:6]
        self.min_voltage_pu = self.ZIP_V[6]

    def get_spu_matrix(self):
        spu_matrix = np.zeros((self.buses.num_elements, 3), dtype=complex)
        for load in self.get_elements():
            bus_idx = self.buses.get_idx(load.related_bus)
            spu_matrix[bus_idx] += load.spu
        return spu_matrix

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

    def _get_zip_val_matrix(self, zip_param=str) -> np.ndarray:
        """
        3 x n matrix of all load.zip_param, aggregated by phase on bus,
        columns indexed by bus
        Zip_params can take any values in
        {'aPQ_p', 'aI_p', 'aZ_p','aPQ_q', 'aI_q', 'aZ_q'}
        """
        param_matrix = np.zeros((self.buses.num_elements, 3))
        for load in self.get_elements():
            load_bus = load.related_bus
            load_ph_idx = load.get_ph_idx_matrix()
            bus_idx = self.buses.get_idx(load_bus)
            param_matrix[bus_idx, load_ph_idx] = getattr(self, zip_param)
        return param_matrix
