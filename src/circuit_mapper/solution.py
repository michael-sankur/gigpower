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
    tolerance = abs((VSLACK[1]) * 10**-9)
    # tolerance = 1e-9
    maxiter = 100
    ZIP_V = [0.10, 0.05, 0.85, 0.10, 0.05, 0.85, 0.80]

    # standardize solution parameter name, index values, columns, and 
    # datatypes across the class
    # see self._init_solution_matrices
    SOLUTION_PARAMS = {
        'V': ['buses', ['A', 'B', 'C'], complex],
        'I': ['lines', ['A', 'B', 'C'], complex],
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
            if element_group == 'lines': # include transformers and vrs
                num_rows = self.circuit.get_total_lines()
            else:
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

    def get_bus_powers(self):
        """
        Total complex powers by bus (load powers and capacitor powers)
        indexed by bus 
        """
        return self.get_load_powers() + self.get_capacitor_powers()

    def calc_Stx(self):
        tx_bus_matrix = self.circuit.get_tx_idx_matrix()
        num_lines = self.circuit.lines.num_elements
        self.Stx = self.V[tx_bus_matrix] * np.conj(self.I)
    
    def calc_Srx(self):
        rx_bus_matrix = self.circuit.get_rx_idx_matrix()
        num_lines = self.circuit.lines.num_elements
        self.Stx = self.V[tx_bus_matrix] * np.conj(self.I)

    def calc_Inode(self) -> None:
        """ Calculate self.Inode (currents consumed at each node) """
        for node in self.network.get_nodes():
            node_V = self.V[node.name]
            node_sV = self.sV[node.name]
            node_I = np.conj(np.divide(node_sV, node_V))
            self.Inode[node.name] = mask_phases(node_I, (3,), node.phases)

    def VMag_df(self):
        """
        returns VMag as a dataframe indexed by node name
        """
        V = self.V_df()
        return V.applymap(lambda cmplx_v: (np.real(cmplx_v)**2 + np.imag(cmplx_v)**2) ** .5)

    def params_df(self):
        """
        returns solution paramaters as a dataframe
        """
        index = ['iterations', 'Vtest', 'Vref', 'tolerance', 'diff']
        data = [self.iterations, self.Vtest,
                self.Vref, self.tolerance, self.diff]
        return pd.DataFrame(data, index).transpose()

    def getLoadPowers(self):
        """
        Return total load powers by bus, calculated from solved V value
        per node.
        """
        data = np.zeros((len(self.network.nodes), 3), dtype=complex)

        for bus_name, bus_idx in self.network.bus_idx_dict.items():
            node = self.network.nodes[bus_name]
            data[bus_idx] = calc_load_power(node, self.V[bus_name])

        return pd.DataFrame(data, self.network.bus_idx_dict.keys(), ['A', 'B', 'C'])

    def getCapPowers(self):
        """
        Return total cap powers by bus, calculated from solved V value
        per node.
        """
        data = np.zeros((len(self.network.nodes), 3), dtype=complex)

        for bus_name, bus_idx in self.network.bus_idx_dict.items():
            node = self.network.nodes[bus_name]
            data[bus_idx] = calc_cap_power(node, self.V[bus_name])

        return pd.DataFrame(data, self.network.bus_idx_dict.keys(), ['A', 'B', 'C'])

    def nomNodePwrs_df(self):
        """
        One time calculation of total nominal node power based on solved V
        equivalent to aP = 1, aI = aQ = 0
        """

        # s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V)).^2) - 1j * cappu + wpu

        data = np.zeros((len(self.network.nodes), 3), dtype=complex)

        for node in self.network.get_nodes():
            node_idx = self.network.bus_idx_dict[node.name]
            nodeV = np.ones((3,), dtype=complex)
            data[node_idx] += calc_total_node_power(
                node, nodeV, [0, 0, 1, 0, 0, 1])
        return pd.DataFrame(data, self.network.bus_idx_dict.keys(), ['A', 'B', 'C'])

    def print_solution(self):
        """
        prints solution to stdout
        """
        print("\n Parameters:")
        print(self.params_df())

        print("\n V solution")
        print(self.V_df())

        print("\n I solution")
        print(self.I_df())

        print("\n Inode solution")
        print(self.Inode_df())

        print("\n Stx solution")
        print(self.Stx_df())

        print("\n Srx solution")
        print(self.Srx_df())

        print("\n sV solution")
        print(self.sV_df())
        print()