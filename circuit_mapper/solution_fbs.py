# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: March 21 2021
# Create FBS Solution class, a namespace for calculations used by fbs

from solution import Solution
import numpy as np


class SolutionFBS(Solution):

    def __init__(self, dss_fp: str):
        super().__init__(dss_fp)  # sets self.circuit
