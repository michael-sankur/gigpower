# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 23 March 2020
# Compare solutions between opendss and SolutiionFBS
# based on tests/test_dss_compare_fbs.py

from circuit_mapper.solution_dss import SolutionDSS
from circuit_mapper.solution_fbs import SolutionFBS
from circuit_mapper.pretty_print import compare_data_frames
from circuit_mapper.circuit import Circuit

import pandas as pd
import numpy as np

import sys
import os
from pathlib import Path
import pytest

DSS_FILE_DIR = Path('./src/nr3_python/')
OUT_DIR = Path('./tests/test_compare_fbs_dss')
OUT_PREFIX = 'FBS_v_DSS_'

GENEROUS = 10e-1
STRICT = 10e-2

Circuit.set_zip_values([0.10, 0.05, 0.85, 0.10, 0.05, 0.85, 0.80])
@pytest.fixture
def dss_solution(dss_file):
    fp = str(Path(DSS_FILE_DIR, dss_file))
    dss_solution = SolutionDSS(str(fp))
    dss_solution.solve()
    return dss_solution

@pytest.fixture
def new_fbs_solution(dss_file):
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
        ('IEEE_13_Bus_allwye_noxfm_noreg.dss', STRICT),  # I fails
        ('IEEE_34_Bus_allwye.dss', GENEROUS),
        ('IEEE_34_Bus_allwye_noxfm_noreg.dss', STRICT),  # not converging
        ('IEEE_37_Bus_allwye.dss', GENEROUS),
        ('IEEE_37_Bus_allwye_noxfm_noreg.dss', STRICT)  # I and sV fail
    ]

)
@pytest.mark.parametrize(
    "solution_param",
    [(param) for param in SolutionFBS.SOLUTION_PARAMS]
)
def test_dss_v_new_fbs(new_fbs_solution, dss_solution, solution_param, 
    dss_file, tolerance):
    """
    Compare the python FBS solution to the opendss solution.
    Writes output to OUTPUT FOLDER.
    """
    fp = Path(DSS_FILE_DIR, dss_file)
    out_file = Path(OUT_DIR, OUT_PREFIX + str(fp.stem) + '_' + solution_param).with_suffix('.out.txt')
    sys.stdout = open(out_file, 'w')
    fbs_vals = new_fbs_solution.get_data_frame(solution_param)
    dss_vals = dss_solution.get_data_frame(solution_param)
    test = ((fbs_vals - dss_vals).abs().max() <= tolerance).all()
    if test:
        print(f"TEST PASSED. FBS CONVERGED = {new_fbs_solution.converged}")
    else:
        print(f"TEST FAILED. FBS CONVERGED = {new_fbs_solution.converged}")
    compare_data_frames(fbs_vals, dss_vals, 'fbs', 'dss', solution_param)
    assert test
