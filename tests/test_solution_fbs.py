# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 21 March 2021
# Test refactor of nr3 module

import numpy as np
import pytest
import pandas as pd

from circuit_mapper.circuit import Circuit
from circuit_mapper.solution_fbs import SolutionFBS

from fbs.fbs.utils import init_from_dss
from fbs.fbs.fbs import fbs, get_solution as get_fbs_solution

DSS_FILE_DIR = 'src/nr3_python/'
DSS_FILE = DSS_FILE_DIR + 'IEEE_13_Bus_allwye.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])

TOLERANCE = 10e-2

@pytest.fixture
def circuit():
    """ map Circuit object once for use in all tests """
    return Circuit(DSS_FILE)


@pytest.fixture
def old_fbs_solution():
    network = init_from_dss(DSS_FILE)
    return get_fbs_solution(fbs(network, 100))


@pytest.fixture
def new_fbs_solution():
    fbs = SolutionFBS(DSS_FILE)
    fbs.maxiter = 100
    fbs.solve()
    return fbs


def test_fbs_V(old_fbs_solution, new_fbs_solution):
    new = new_fbs_solution.get_data_frame('V')
    old = old_fbs_solution.V_df()
    print(new)
    print(old)
    assert ((new - old).abs().max() <= TOLERANCE).all()
    # assert new.equals(old)
    
def test_fbs_I(old_fbs_solution, new_fbs_solution):
    new = new_fbs_solution.get_data_frame('I')
    old = old_fbs_solution.I_df()
    print(new)
    print(old)
    assert ((new - old).abs().max() <= TOLERANCE).all()
    # assert new.equals(old)


def test_fbs_sV(old_fbs_solution, new_fbs_solution):
    new = new_fbs_solution.get_data_frame('sV')
    old = old_fbs_solution.sV_df()
    print(new)
    print(old)
    assert ((new - old).abs().max() <= TOLERANCE).all()
    # assert new.equals(old)