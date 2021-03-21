# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 21 March 2021
# Test refactor of nr3 module

import numpy as np
import pytest

from circuit import Circuit
from solution_fbs import SolutionFBS

from fbs.utils import init_from_dss
from fbs.fbs import fbs, get_solution as get_fbs_solution

LOCAL_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
DSS_FILE = LOCAL_DIR + 'IEEE_13_Node/IEEE_13_Bus_allwye.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])


@pytest.fixture
def circuit():
    """ map Circuit object once for use in all tests """
    return Circuit(DSS_FILE)


@pytest.fixture
def old_fbs_solution():
    network = init_from_dss(DSS_FILE)
    return get_fbs_solution(fbs(network))


@pytest.fixture
def new_fbs_solution():
    return SolutionFBS(DSS_FILE)


def test_fbs_V(old_fbs_solution, new_fbs_solution):
    assert (old_fbs_solution.V_df() == new_fbs_solution.V_df()).all()


def test_fbs_I(old_fbs_solution, new_fbs_solution):
    assert (old_fbs_solution.I_df() == new_fbs_solution.I_df()).all()


def test_fbs_Stx(old_fbs_solution, new_fbs_solution):
    assert (old_fbs_solution.Stx_df() == new_fbs_solution.Stx_df()).all()


def test_fbs_Srx(old_fbs_solution, new_fbs_solution):
    assert (old_fbs_solution.Srx_df() == new_fbs_solution.Srx_df()).all()


def test_fbs_sV(old_fbs_solution, new_fbs_solution):
    assert (old_fbs_solution.sV_df() == new_fbs_solution.sV()).all()
