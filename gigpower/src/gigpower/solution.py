# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Solution superclass

from . circuit import Circuit
from . volt_var_controller import VoltVARController
from typing import Iterable, Dict, Any
import opendssdirect as dss

import re
import numpy as np
import pandas as pd


class Solution():

    # class variables set for all SolutionNR3 instances
    # TODO: If any of these need to be set by instance, move into self.__init__
    SLACKIDX = 0  # assume slack bus is at index 0

    # TODO: VSLACK the same for all objects. Write a SETVSLACK method on the class.
    VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)],
                      dtype=complex)

    maxiter = 100

    # standardize solution parameter name, index values, columns, and
    # datatypes across the class
    # see self._init_solution_matrices
    SOLUTION_PARAMS = {
        'V': ['buses', ['A', 'B', 'C'], complex],
        'I': ['lines', ['A', 'B', 'C'], complex],
        'sV': ['buses', ['A', 'B', 'C'], complex],
        'Vmag': ['buses', ['A', 'B', 'C'], float],
        'Stx': ['lines', ['A', 'B', 'C'], complex],
        'Srx': ['lines', ['A', 'B', 'C'], complex]
    }

    @classmethod
    def set_zip_values(cls, zip_v):
        """
        sets zip values for the Solution class
        param zip_V: List or nd.array with 7 values
        [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min voltage pu]
        Note that zip values are set both on the Solution class and Circuit
        class
        """
        cls.ZIP_V = np.asarray(zip_v)
        cls.aZ_p, cls.aI_p, cls.aPQ_p = cls.ZIP_V[0:3]
        cls.aZ_q, cls.aI_q, cls.aPQ_q = cls.ZIP_V[3:6]
        cls.min_voltage_pu = cls.ZIP_V[6]
        Circuit._set_zip_values(zip_v)

    # TODO: Make a 'solution.set_tolerance()' method
    def __init__(self, dss_fp: str, zip_v: np.ndarray = np.asarray([
                                                                    0.10, 0.05,
                                                                    0.85, 0.10,
                                                                    0.05, 0.85,
                                                                    0.80])):
        """
        sets up a Solution object with a pointer to a Circuit mapped from opendss
        Solutions keep a pointer to the dss object used to map the Circuit
        param dss_fp: path to the dss file
        param zip_V: optional Zip Values to set for all Solutions and Circuits
        defaults to [.1, .05, .85, .1, .05. ,.85, .8]
        """
        self.set_zip_values(zip_v)
        #  setup calls to opendss----------------------------------------------
        self.dss = dss  # save dss instance to object
        dss.run_command('Redirect ' + dss_fp)
        dss.Solution.Solve()  # solve first for base values
        # set zip values for all Solutions and all Circuits

        dss.Solution.Solve()  # solve again to set zip values on dss

        #  map Circuit -----------------------------------------
        self.circuit = Circuit(dss)

        # default orientation, n x 3 or 3 x n
        # consider making this a class variable
        self._orient = 'rows'
        self.circuit._orient = self._orient

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
        # cache sV params for load calcualtions
        self._set_sV_params()

    def _set_sV_params(self):
        self.phase_matrix = self.circuit.buses.get_phase_matrix(self._orient)
        self.spu = self.circuit.get_spu_matrix()
        self.aPQ = self.circuit.get_aPQ_matrix()
        self.aI = self.circuit.get_aI_matrix()
        self.aZ = self.circuit.get_aZ_matrix()
        self.cappu = self.circuit.get_cappu_matrix()
        self.wpu = self.circuit.get_wpu_matrix()

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
            # include transformers and vrs, for fbs
            if element_group == 'lines':
                num_rows = self.circuit.get_total_lines()
            else:
                num_rows = getattr(getattr(self.circuit, element_group), 'num_elements')
            setattr(self, param, np.zeros((num_rows, len(cols)), dtype=datatype))

    def get_data_frame(self, param: str, orient: str = '') -> pd.DataFrame:
        """
        Returns a DataFrame for the specified solution paramater.
        param: must be in SOLUTION_PARAMS
        orient: optional, 'rows' or 'cols', defaults to self._orient
        """
        if not orient:
            orient = self._orient
        try:
            element_group, cols, data_type = self.__class__.SOLUTION_PARAMS.get(param)
            # force a deep copy to avoid pointer issues
            index = [_ for _ in getattr(self.circuit, element_group).all_names()]
            data = getattr(self, param)
            # FBS saves transformers and voltage regulators in self.I
            # ignore transformers and voltage regulators for solution values
            if element_group == 'lines' and self.__class__.__name__ == 'SolutionFBS':
                data = data[0:self.circuit.lines.num_elements]
            if self.__class__.__name__ == 'SolutionNR3':
                data = data.transpose()            
            if orient == 'cols':
                data = data.transpose()
                # force a deep copy swap to avoid pointer issues
                temp = [_ for _ in cols]
                cols = [_ for _ in index]
                index = temp
            return pd.DataFrame(data=data, index=index, columns=cols, dtype=data_type)
        except KeyError:
            print(f"Not a valid solution parameter. Valid parameters: \
                  {self.__class__.SOLUTION_PARAMS.keys()}")

    def get_nominal_bus_powers(self, orient: str = ''):
        """
        Returns a DataFrame for self.Circuit's powers by bus
        param: must be in SOLUTION_PARAMS
        orient: optional, 'rows' or 'cols', defaults to self._orient
        """
        data = self.circuit.get_nominal_bus_powers() * 1000
        index = self.circuit.buses.all_names()
        cols = (['A', 'B', 'C'])
        data_type = complex
        if self.__class__.__name__ == 'SolutionNR3':
            data = data.transpose()
        if orient == 'cols':
            data = data.transpose()
            # force a deep copy swap to avoid pointer issues
            temp = [_ for _ in cols]
            cols = [_ for _ in index]
            index = temp
        return pd.DataFrame(data=data, index=index, columns=cols, dtype=data_type)

    def calc_sV(self, bus=None):
        """
        Used to calculate total node powers, accounting for loads
        and capacitors, oriented according to self._orient
        self._orient = 'rows' -> returns n x 3
        self._orient = 'cols' -> returns 3 x n
        param bus: if a Bus is given, updates self.sV only at the index of the Bus
        """
        V = self.V
        spu = self.spu
        aPQ, aI, aZ = self.aPQ, self.aI, self.aZ
        cappu, wpu = self.cappu, self.wpu
        phase_matrix = self.phase_matrix

        if bus:
            bus_idx = self.circuit.buses.get_idx(bus)
            V = self.V[bus_idx]
            spu = self.spu[bus_idx]
            aPQ, aI, aZ = self.aPQ[bus_idx], self.aI[bus_idx], self.aZ[bus_idx]
            cappu, wpu = self.cappu[bus_idx], self.wpu[bus_idx]
            phase_matrix = self.phase_matrix[bus_idx]

        # TODO: confirm if cappu.real needs to be multiplied by abs(V)**2
        # nr3 map_solution does not do that, but fbs requires it for
        # results to be consistent with opendss
        update = spu * (aPQ + aI * np.abs(V) + aZ * np.abs(V) ** 2) - \
            1j * cappu * np.abs(V)**2 + 1j * wpu

        update[phase_matrix == 0] = 0

        if bus:
            self.sV[bus_idx] = update
        else:
            if self._orient == 'rows':
                outer, inner = self.circuit.buses.num_elements, 3
            elif self._orient == 'cols':
                outer, inner = 3, self.circuit.buses.num_elements
            # for i in range(outer):
            #     for j in range(inner):
            #         if np.abs(update[i, j].real) <= 1e-12:
            #             update[i, j] = 0 + update[i, j].imag
            #         if np.abs(self.sV[i, j].imag) <= 1e-12:
            #             update[i, j] = update[i, j].real + 0
            self.sV = update

    def calc_Stx(self):
        num_lines = self.circuit.lines.num_elements
        buses = self.circuit.get_tx_idx_matrix()[0:num_lines]
        self.Stx = self.V[buses] * np.conj(self.I[0:num_lines])

    def calc_Srx(self):
        num_lines = self.circuit.lines.num_elements
        buses = self.circuit.get_rx_idx_matrix()[0:num_lines]
        self.Srx = self.V[buses] * np.conj(self.I[0:num_lines])

    def calc_Vmag(self) -> np.ndarray:
        self.Vmag = abs(self.V)

    @classmethod
    def get_params(cls) -> Dict:
        """
        returns solution paramaters as a dictionary
        """
        params = {
            'slack_idx': cls.SLACKIDX,
            'v_slack': cls.VSLACK,
            'max_iter': cls.maxiter,
            'zip_v': cls.ZIP_V
        }
        return params

    def _set_orient(self, orient: str):
        """
        Sets matrix orientation on self and self.circuit, indicating whether 
        shapes are n x 3 (rows) or 3 x n (columns)
        param orient: 'rows' or 'cols'
        TODO: implement this to change all matrices currently set on solution 
        to the correct orientation
        """
        self._orient = orient
        self.circuit._orient = orient
        pass

    def get_V(self, orient='') -> pd.DataFrame:
        """ Return solved voltage as a DataFrame"""
        return self.get_data_frame('V', orient)

    def get_Vmag(self, orient='') -> pd.DataFrame:
        """ Return solved voltage magnitude as a DataFrame"""
        return self.get_data_frame('Vmag', orient)

    def get_I(self, orient='') -> pd.DataFrame:
        """ Return solved line currents as a DataFrame"""
        return self.get_data_frame('I', orient)

    def get_Stx(self, orient='') -> pd.DataFrame:
        """ Return solved incoming line powers as a DataFrame"""
        return self.get_data_frame('Stx', orient)

    def get_Srx(self, orient='') -> pd.DataFrame:
        """ Return solved outgoing line powers as a DataFrame"""
        return self.get_data_frame('Srx', orient)

    def get_sV(self, orient='') -> pd.DataFrame:
        """ Return solved total bus powers as a DataFrame"""
        return self.get_data_frame('sV', orient)

    def print_solution(self):
        """ Prints solution values to stdout, row-major order."""
        pd.options.display.float_format = '{:,.3f}'.format

        print("\nParameters:")
        for p, v in self.get_params().items():
            print(f'{p}: {v}')
        print('-' * 80)

        print("\nV solution")
        print(self.get_V('rows'))
        print('-' * 80)

        print("\nVmag solution")
        print(self.get_Vmag('rows'))
        print('-' * 80)

        print("\nI solution")
        print(self.get_I('rows'))
        print('-' * 80)

        print("\nStx solution")
        print(self.get_Stx('rows'))
        print('-' * 80)

        print("\nSrx solution")
        print(self.get_Srx('rows'))
        print('-' * 80)

        print("\nsV solution")
        print(self.get_sV('rows'))
        print('-' * 80)
