# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Solution superclass

from circuit import Circuit
import utils
import opendssdirect as dss
ZIPV = [0.10, 0.05, 0.85, 0.10, 0.05, 0.85, 0.80]


class Solution():

    def __init__(self, dss_fp: str):
        """
        sets up a solution object with circuit mapped from opendss
        """
        self.dss = dss.run_command('Redirect ' + dss_fp)
        dss.Solution.Solve()  # solve first for base values
        utils.set_zip_values(dss, ZIPV)  # set zip values before mapping Circuit!
        dss.Solution.Solve()  # solve again to set zip values
        self.circuit = Circuit(dss)
