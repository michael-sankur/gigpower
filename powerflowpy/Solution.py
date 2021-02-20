# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Solution superclass

import Circuit


class Solution():

    def __init__(self, dss) -> Solution:
        """
        sets up a solution object given an opendss object's current state.
        make sure that zip values are set up
        """
        self.circuit = Circuit(dss)
