# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 23 March 2020
# Compare solutions between opendss and SolutiionFBS
# based on tests/test_dss_compare_fbs.py

from circuit_mapper.solution_dss import SolutionDSS
from circuit_mapper.solution_fbs import SolutionFBS
from circuit_mapper.solution_nr3 import SolutionNR3
from circuit_mapper.pretty_print import compare_data_frames
from circuit_mapper.circuit import Circuit

import numpy as np

import sys
import os
from pathlib import Path
import pytest

DSS_FILE_DIR = Path('./src/nr3_python/')
OUT_DIR_FBS = Path('./tests/test_compare_fbs_dss')
OUT_DIR_NR3 = Path('./tests/test_compare_nr3_dss')

# thresholds for passing tests
GENEROUS = .5
STRICT = .01


@pytest.fixture
def dss_solution(dss_file):
    fp = str(Path(DSS_FILE_DIR, dss_file))
    dss_solution = SolutionDSS(str(fp))
    dss_solution.solve()
    return dss_solution


def new_fbs_solution(dss_file):
    fp = str(Path(DSS_FILE_DIR, dss_file))
    solution = SolutionFBS(fp)
    solution.solve()
    return solution


def new_nr3_solution(dss_file):
    fp = str(Path(DSS_FILE_DIR, dss_file))
    solution = SolutionNR3(fp)
    solution.solve()
    return solution


def setup_module():
    if not OUT_DIR_NR3.exists():
        os.makedirs(OUT_DIR_NR3)
    if not OUT_DIR_FBS.exists():
        os.makedirs(OUT_DIR_FBS)


@pytest.mark.parametrize(
    "algorithm,out_dir",
    [
        ("NR3", OUT_DIR_NR3),
        ("FBS", OUT_DIR_FBS),
    ]
)
@pytest.mark.parametrize(
    "zip_values,zip_name",
    [
        (np.asarray([1, 0, 0, 1, 0, 0, .8]),  'constant_impedance'),
        (np.asarray([0, 0, 1, 0, 0, 1, .8]),  'constant_power'),
        (np.asarray([0.10, 0.05, 0.85, 0.10, 0.05,
                     0.85, 0.80]), 'test_zip')
    ]
)
@pytest.mark.parametrize(
    "dss_file,tolerance",
    [
        ('IEEE_13_Bus_allwye.dss', GENEROUS),
        ('IEEE_13_Bus_allwye_noxfm_noreg.dss', STRICT),
        ('IEEE_34_Bus_allwye.dss', GENEROUS),
        ('IEEE_34_Bus_allwye_noxfm_noreg.dss', STRICT),
        ('IEEE_37_Bus_allwye.dss', GENEROUS),
        ('IEEE_37_Bus_allwye_noxfm_noreg.dss', STRICT)
    ]

)
@pytest.mark.parametrize(
    "solution_param",
    [(param) for param in SolutionFBS.SOLUTION_PARAMS]
)
def test_dss_v_new_fbs(algorithm, out_dir, dss_solution, zip_values,
                       zip_name, solution_param, dss_file, tolerance):
    """
    Compare the python FBS solution to the opendss solution.
    Writes output to OUTPUT FOLDER.
    """
    if algorithm == 'NR3':
        new_solution = new_nr3_solution(dss_file)
    elif algorithm == 'FBS':
        new_solution = new_fbs_solution(dss_file)
    Circuit.set_zip_values(zip_values)
    fp = Path(DSS_FILE_DIR, dss_file)
    out_prefix = f"{algorithm}_v_DSS_"
    out_file = Path(out_dir, out_prefix + str(fp.stem) + '_' +
                    zip_name + '-' +
                    solution_param).with_suffix('.out.txt')
    sys.stdout = open(out_file, 'w')
    new_vals = new_solution.get_data_frame(solution_param)
    dss_vals = dss_solution.get_data_frame(solution_param)

    # set sV pass threshold to 5% of opendss sV magnitude
    if solution_param == 'sV':
        dss_thresholds = dss_vals.abs() * .05
        diffs = (new_vals.abs() - dss_vals.abs()).abs()
        test = ((diffs <= dss_thresholds).max()).all()
    else:
        if solution_param == 'Vmag': # set Vmag tolerance to strict in all cases
            tolerance = STRICT
        test = ((new_vals - dss_vals).abs().max() <= tolerance).all()
    if test:
        print(f"TEST PASSED. {algorithm} CONVERGED = {new_solution.converged}")
    else:
        print(f"TEST FAILED. {algorithm} CONVERGED = {new_solution.converged}")
    print(f"Zip values = {zip_values}")
    compare_data_frames(new_vals, dss_vals, algorithm, 'dss', solution_param)
    assert test
