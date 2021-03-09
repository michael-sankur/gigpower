# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021

import numpy as np  # type: ignore
from typing import List
from load import Load
from capacitor import Capacitor
from circuit_element import CircuitElement


class Bus(CircuitElement):
    dss_module_name = 'Bus'

    def __init__(self, name: str, dss) -> None:
        super().__init__(name, dss)
        self.upstream_buses = []
        self.downstream_buses = []
        self.loads: List[Load] = []
        self.capacitors: List[Capacitor] = []
        self.sum_spu = np.zeros(3)  # 3x1 complex array that holds the sum of all load.spu on this node
        self.sum_cappu = np.zeros(3)  # 3x1 array that holds the sum of all capacitor.cappu on this node

    def _set_related_bus(self, dss):
        self.related_bus = self.__name__

    # def _set_phases_from_bus(self, dss, bus_name):
    #     """ override superclass to set phases  """
    #     self.phases = []

    def _set_spu(self, subtract_spu, add_spu):
        """set the bus' sum spu by subtracting subtract_spu, and adding add_spu"""
        self.sum_spu = -1 * subtract_spu + add_spu
