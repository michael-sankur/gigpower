# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: March 21 2021
# Create FBS Solution class, a namespace for calculations used by fbs

from solution import Solution
import numpy as np


class SolutionFBS(Solution):

    def __init__(self, dss_fp: str):
        super().__init__(dss_fp)  # sets self.circuit
        self._init_XNR()
        self._init_slack_bus_matrices()
        self._init_KVL_matrices()
        self._init_KCL_matrices()
        self._init_KVL_matrices_vregs()
