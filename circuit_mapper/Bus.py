# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021

import numpy as np  # type: ignore
from typing import List
import Load
import Capacitor


class Bus(CircuitElement):

    def __init__(self, name: str, dss) -> None:
        super().__init__(self, name)
        self.upstream_buses = []
        self.downstream_buses = []
        self.loads: List[Load] = []
        self.capacitors: List[Capacitor] = []
        self.sum_spu = np.zeros(3)  # 3x1 complex array that holds the sum of all load.spu on this node
        self.sum_cappu = np.zeros(3)  # 3x1 array that holds the sum of all capacitor.cappu on this node
        self.Vbase = 0.0
        self.Ibase = 0.0
        self.Zbase = 0.0
