# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 21 March 2021
# Test refactor of nr3 module

import numpy as np
import pytest
import pandas as pd

from circuit_mapper.circuit import Circuit
from circuit_mapper.solution_fbs import SolutionFBS
from circuit_mapper.pretty_print import compare_data_frames

from fbs.fbs.utils import init_from_dss
from fbs.fbs.fbs import fbs, get_solution as get_fbs_solution
import sys
import os
from pathlib import Path

DSS_FILE_DIR = Path('./src/nr3_python/')
OUT_DIR = Path('./tests/test_compare_old_fbs')
OUT_PREFIX = 'FBS_v_FBS_'

GENEROUS = 5e-1
STRICT = 1e-1

@pytest.fixture
def old_fbs_solution(dss_file):
    # dss_file = request.param
    fp = str(Path(DSS_FILE_DIR, dss_file))
    network = init_from_dss(fp)
    return get_fbs_solution(fbs(network))

@pytest.fixture
def new_fbs_solution(dss_file):
    # dss_file = request.param
    fp = str(Path(DSS_FILE_DIR, dss_file))
    solution = SolutionFBS(fp)
    solution.solve()
    return solution

def setup_module():
    if not OUT_DIR.exists():
        os.makedirs(OUT_DIR)

@pytest.mark.parametrize(
    "dss_file,tolerance",
    [
        ('IEEE_13_Bus_allwye.dss', GENEROUS),
        ('IEEE_13_Bus_allwye_noxfm_noreg.dss', GENEROUS),
        ('IEEE_34_Bus_allwye.dss', GENEROUS),
        ('IEEE_34_Bus_allwye_noxfm_noreg.dss', GENEROUS),
        ('IEEE_37_Bus_allwye.dss', GENEROUS),
        ('IEEE_37_Bus_allwye_noxfm_noreg.dss', GENEROUS)
    ],
)
@pytest.mark.parametrize(
    "new_fbs_param, old_fbs_method",
    [
       ('V', 'V_df'),
       ('I', 'I_df'),
       ('sV', 'sV_df'),
    ],
)

def test_fbs_v_fbs(new_fbs_solution, old_fbs_solution, dss_file, tolerance, new_fbs_param,
    old_fbs_method):
    """ Test all Files X (V, I, sV). Save log files to OUT_DIR"""
    fp = Path(DSS_FILE_DIR, dss_file)
    new = new_fbs_solution.get_data_frame(new_fbs_param)
    old = getattr(old_fbs_solution, old_fbs_method)()
    out_file = Path(OUT_DIR, new_fbs_param + '_' + OUT_PREFIX + str(fp.stem)).with_suffix('.out.txt')
    sys.stdout = open(out_file, 'w')
    test = ((new - old).abs().max() <= tolerance).all()
    if test:
        print("TEST PASSED")
    else:
        print("TEST FAILED")
    compare_data_frames(new, old, 'new_fbs', 'old_fbs', new_fbs_param)
    assert test

