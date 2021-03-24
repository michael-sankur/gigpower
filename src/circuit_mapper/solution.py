# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Solution superclass

from . circuit import Circuit
from . volt_var_controller import VoltVARController
from . utils import set_zip_values
import opendssdirect as dss

import re
import numpy as np
from typing import Dict
import pandas as pd


class Solution():

    # class variables set for all SolutionNR3 instances
    # TODO: If any of these need to be set by instance, move into self.__init__
    SLACKIDX = 0  # assume slack bus is at index 0
    VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])
    V0, I0 = None, None

    # # TODO: Make a better 'solution.set_tolerance(ref_node, error)' method
    # set tolerance with phase B reference voltage
    # tolerance = abs((solution.Vref[1]) * 10**-9)
    tolerance = 1e-9
    maxiter = 100
    ZIP_V = [0.10, 0.05, 0.85, 0.10, 0.05, 0.85, 0.80]

    # standardize solution parameter name, index values, columns, and 
    # datatypes across the class
    # see self._init_solution_matrices
    SOLUTION_PARAMS = {
        'V': ['buses', ['A', 'B', 'C'], complex],
        'I': ['lines', ['A', 'B', 'C'], complex],
        'Inode': ['buses', ['A', 'B', 'C'], complex],
        'Stx': ['lines', ['A', 'B', 'C'], complex],
        'Srx': ['lines', ['A', 'B', 'C'],  complex],
        'sV': ['buses', ['A', 'B', 'C'], complex]}

    def __init__(self, dss_fp: str):
        """
        sets up a Solution object with a pointer to a Circuit mapped from opendss
        Solutions keep a pointer to the dss object used to map the Circuit
        """
        #  setup calls to opendss----------------------------------------------
        self.dss = dss.run_command('Redirect ' + dss_fp)
        dss.Solution.Solve()  # solve first for base values
        set_zip_values(dss, self.__class__.ZIP_V)
        dss.Solution.Solve()  # solve again to set zip values

        #  map Circuit and vvc objects-----------------------------------------
        self.circuit = Circuit(dss)
        self.volt_var_controllers = self.parse_vvc_objects(dss_fp)

        #  initialize solution parameters---------------------------------------
        self.iterations = 0
        # stores the tolerance at most recent completed iteration
        self.solution_tolerance = -1
        # stores the final value of Vtest - Vref at convergence
        self.convergence_diff = -1  
        self._init_solution_matrices()

        # Voltage parameters. TODO: are these only for fbs?
        # If so, move to solution_fbs.py
        self.Vtest = np.zeros(3, dtype='complex')
        self.Vref = np.array(
            [1, np.exp(1j*240*np.pi/180), np.exp(1j*120*np.pi/180)], dtype=complex)

    def _init_solution_matrices(self):
        """
        Initializes matrices to store solution values in ndarrays as follows:
        V: num_buses x 3, complex pu voltage, by bus index
        I: num_lines x 3, complex pu current phasors, by line index
        Inode: num_buses x 3, complex pu current phasors delivered to bus, by bus index
        Stx: num_lines x 3, line transmitting end power, by line index
        Srx: num_lines x 3, line receiving end power, by line index
        sV: num_buses x 3, total powers at each bus, by bus index
        """
        for param in self.__class__.SOLUTION_PARAMS:
            element_group, cols, datatype = self.__class__.SOLUTION_PARAMS[param]
            num_rows = getattr(getattr(self.circuit, element_group), 'num_elements')
            setattr(self, param, np.zeros((num_rows, len(cols)), dtype=datatype))

    def get_data_frame(self, param: str) -> pd.DataFrame:
        """
        Returns a DataFrame for the specified solution paramater.
        param: must be in SOLUTION_PARAMS
        """
        try:
            element_group, cols, data_type = self.__class__.SOLUTION_PARAMS.get(param)
            index = getattr(self.circuit, element_group).all_names()
            data = getattr(self, param)
            return pd.DataFrame(data=data, index=index, columns=cols, dtype=data_type)
        except KeyError:
            print(f"Not a valid solution parameter. Valid parameters: \
                  {self.__class__.SOLUTION_PARAMS.keys()}")

    def parse_vvc_objects(self, fn: str):
        """ From 20180601/PYTHON/lib/dss_vvc.py by @kathleenchang"""
        # Parse VVC lines in DSS file
        vvarobjects = []
        f = open(fn, "r")

        for l in f:
            if re.findall('(?i)New VVC', l):
                bp = []

                #array of VVC breakpoints
                breakpoints = str.split(re.findall(r"(?i)BP=\s*([^\n\r]*)", l)[0], ',')

                for elem in breakpoints:
                    point = re.findall("[0-9.]*",  elem)
                    for i in point:
                        if i:
                            bp.append(float(i))
                print(bp)

                # zero-indexed phase, is one-indexed in DSS file
                phase = int(re.findall('(?i)phase=\s*([0-9]*)', l)[0]) - 1
                print(phase)

                minkvar = float(re.findall('(?i)min_kvar=([-0-9]*)', l)[0])
                print(minkvar)

                maxkvar = float(re.findall('(?i)max_kvar=([-0-9]*)', l)[0])
                print(maxkvar)

                bus = re.findall(r"(?i)bus=([\w.]*)\s", l)[0]
                bus = re.findall(r"[\w]*", bus)[0]
                print(bus)

                # create volt var object
                voltvarobject = VoltVARController(bp, minkvar, maxkvar, bus, phase)
                vvarobjects.append(voltvarobject)
                print("\n --------------")

        for e in vvarobjects:
            print(e)
        return vvarobjects
