# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 21 March 2021
# Test refactor of nr3 module

import numpy as np
import opendssdirect as dss
import pytest

from circuit import Circuit
from solution_nr3 import SolutionNR3

# current nr3 dependencies
import sys
sys.path.append('/Users/elainelaguerta/Dropbox/LBNL/LinDist3Flow/20180601/PYTHON')
from lib.DSS_parameters import relevant_openDSS_parameters
from lib.basematrices import basematrices
from lib.helper import transformer_regulator_parameters, nominal_load_values, cap_arr

LOCAL_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
DSS_FILE = LOCAL_DIR + 'IEEE_13_Node/IEEE_13_Bus_allwye.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])


@pytest.fixture
def circuit():
    """ map Circuit object once for use in all tests """
    return Circuit(DSS_FILE)


@pytest.fixture
def nr3_solution():
    return SolutionNR3(DSS_FILE)
