# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 21 March 2021
# Test refactor of nr3 module

import numpy as np
import opendssdirect as dss
import pytest

from gigpower.circuit import Circuit
from gigpower.solution_nr3 import SolutionNR3

# current nr3 dependencies
from nr3_python.lib.DSS_parameters import relevant_openDSS_parameters
from nr3_python.lib.basematrices import basematrices
from nr3_python.lib.helper import transformer_regulator_parameters, nominal_load_values, cap_arr
from nr3_python.lib.NR3 import NR3

LOCAL_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
DSS_FILE = LOCAL_DIR + 'IEEE_13_Node/IEEE_13_Bus_allwye.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])


@pytest.fixture
def new_nr3_object():
    return SolutionNR3(DSS_FILE)


@pytest.fixture
def old_nr3_object():
    pass
    # TODO: create an nr3 solution from lib.NR3

# NR3 VVC TESTS---------------------------------------------------------------------


@pytest.mark.xfail
def test_get_Q(circuit, circuit_object, nr3_object, voltage):
    nr_Q = nr3_object.get_Q(voltage)
    ckt_Q = circuit_object.get_Q(voltage)
    assert(nr_Q == ckt_Q)


@pytest.mark.xfail
def test_get_prevQ(circuit, circuit_object, nr3_object):
    nr_prevQ = nr3_object.get_prevQ()
    ckt_prevQ = circuit_object.get_prevQ()

    assert(nr_prevQ == ckt_prevQ)


@pytest.mark.xfail
def test_prevQ_list(circuit, circuit_object, nr3_object):
    nr_list = nr3_object.get_prevQ_list()
    ckt_list = circuit_object.get_prevQ_list()
    assert((nr_list == ckt_list).all)


@pytest.mark.xfail
def test_vvc_parameters(circuit, circuit_object, nr3_object):
    nr_busName = nr3_object.get_busName()
    nr_phase = nr3_object.get_phase()

    ckt_busName = circuit.get_vvc_busname(circuit_object)
    ckt_phase = circuit.get_vvc_phase(circuit_object)
    assert(nr_busName == ckt_busName)
    assert(nr_phase == ckt_phase)
