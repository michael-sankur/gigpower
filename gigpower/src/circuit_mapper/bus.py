# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021

import numpy as np  # type: ignore
from typing import List
from . load import Load
from . capacitor import Capacitor
from . circuit_element import CircuitElement


class Bus(CircuitElement):
    dss_module_name = 'Bus'

    def __init__(self, name: str, dss) -> None:
        super().__init__(name, dss)
        self.upstream_buses = []
        self.downstream_buses = []
        self.loads: List[Load] = []
        self.capacitors: List[Capacitor] = []

    def _set_related_bus(self, dss):
        self.related_bus = self.__name__

    def _set_spu(self, subtract_spu, add_spu):
        """set the bus' sum spu by subtracting subtract_spu, and adding add_spu"""
        self.sum_spu = -1 * subtract_spu + add_spu
