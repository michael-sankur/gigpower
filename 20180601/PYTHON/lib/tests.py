

import numpy as np
import opendssdirect as dss
import pytest

from circuit import Circuit
from solution_nr3 import SolutionNR3
from solution import Solution

# current nr3 dependencies
import sys
sys.path.append('/Users/elainelaguerta/Dropbox/LBNL/LinDist3Flow/20180601/PYTHON')
from lib.vvc import voltVARControl #???

LOCAL_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
DSS_FILE = LOCAL_DIR + 'IEEE_13_Node/IEEE_13_Bus_allwye.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])


@pytest.fixture
def circuit():
    """ map Circuit object once for use in all tests """
    dss.run_command('Redirect ' + DSS_FILE)
    return Circuit(dss)


# @pytest.fixture
# def nr3_solution():
#     return SolutionNR3(DSS_FILE)


# @pytest.fixture
# def nr3_DSS_parameters():
#     return relevant_openDSS_parameters(DSS_FILE, -1)


# @pytest.fixture
# def xfm_vr_parameters():
#     return transformer_regulator_parameters()


# @pytest.fixture
# def nr3_basematrices():
#     slack_idx = 0
#     Vslack = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])
#     V0, I0 = None, None
#     return basematrices(DSS_FILE, slack_idx, Vslack, V0, I0)


# NR3 VVC TESTS---------------------------------------------------------------------
def test_get_Q(circuit, circuit_object, nr3_object, voltage):
    nr_Q = nr3_object.get_Q(voltage)
    ckt_Q = circuit_object.get_Q(voltage)
    assert(nr_Q == ckt_Q)

def test_get_prevQ(circuit, circuit_object, nr3_object):
    nr_prevQ = nr3_object.get_prevQ()
    ckt_prevQ = circuit_object.get_prevQ()

    assert(nr_prevQ == ckt_prevQ)

def test_prevQ_list(circuit, circuit_object, nr3_object):
    nr_list = nr3_object.get_prevQ_list()
    ckt_list = circuit_object.get_prevQ_list()
    assert((nr_list == ckt_list).all)

def test_vvc_parameters(circuit, circuit_object, nr3_object):
    nr_busName = nr3_object.get_busName()
    nr_phase = nr3_object.get_phase()

    ckt_busName = circuit.get_vvc_busname(circuit_object)
    ckt_phase = circuit.get_vvc_phase(circuit_object)
    assert(nr_busName == ckt_busName)
    assert(nr_phase == ckt_phase)

    