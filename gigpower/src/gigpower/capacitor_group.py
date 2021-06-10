from gigpower.circuit_element_group import CircuitElementGroup
from . capacitor import Capacitor
import numpy as np
import pandas as pd

class CapacitorGroup(CircuitElementGroup):
    dss_module_name = 'Capacitors'
    ele_class = Capacitor
 
    def __init__(self, dss, bus_group):
        self.buses = bus_group
        super().__init__(dss)

    def get_cappu_matrix(self):
        """
        return 3 x n matrix of kvar values summed over Bus
        columns indexed by Bus index, padded by phase
        """
        if self.num_elements > 0:
            return self._get_attr_by_bus('cappu', orient='cols')
        else:
            return np.zeros((self.buses.num_elements, 3), dtype=float)
