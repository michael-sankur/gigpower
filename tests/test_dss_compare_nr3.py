# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 29 March 2020
# Compare solutions between opendss and Solution NR3
# based on tests/test_dss_compare_fbs.py

from circuit_mapper.solution_dss import SolutionDSS
from circuit_mapper.solution_nr3 import SolutionNR3
from circuit_mapper.circuit import Circuit
from circuit_mapper.pretty_print import compare_data_frames

import sys
import os
from pathlib import Path
import pytest

DSS_FILE_DIR = Path('./src/nr3_python/')
OUT_DIR = Path('./tests/test_compare_nr3_dss')
OUT_PREFIX = 'NR3_v_DSS_'

GENEROUS = 10e-1
STRICT = 10e-2


@pytest.fixture
def dss_solution(dss_file):
    fp = str(Path(DSS_FILE_DIR, dss_file))
    dss_solution = SolutionDSS(str(fp))
    dss_solution.circuit.set_zip_values([1, 0, 0, 1, 0, 0, 0.80])
    dss_solution.solve()
    return dss_solution


@pytest.fixture
def new_nr3_solution(dss_file):
    fp = str(Path(DSS_FILE_DIR, dss_file))
    solution = SolutionNR3(fp)
    solution.circuit.set_zip_values([1, 0, 0, 1, 0, 0, 0.80])
    solution.solve()
    return solution


def setup_module():
    if not OUT_DIR.exists():
        os.makedirs(OUT_DIR)


@pytest.mark.parametrize(
    "dss_file,tolerance",
    [
        # ('IEEE_13_Bus_allwye.dss', GENEROUS),
        ('IEEE_13_Bus_allwye_noxfm_noreg.dss', STRICT),
        # ('IEEE_34_Bus_allwye.dss', GENEROUS),
        # ('IEEE_34_Bus_allwye_noxfm_noreg.dss', STRICT),
        # ('IEEE_37_Bus_allwye.dss', GENEROUS),
        # ('IEEE_37_Bus_allwye_noxfm_noreg.dss', STRICT)
    ]

)
@pytest.mark.parametrize(
    "solution_param",
    [(param) for param in SolutionNR3.SOLUTION_PARAMS]
)
def test_dss_v_new_nr3(new_nr3_solution, dss_solution, solution_param,
                       dss_file, tolerance):
    """
    Compare the python FBS solution to the opendss solution.
    Writes output to OUTPUT FOLDER.
    """
    fp = Path(DSS_FILE_DIR, dss_file)
    out_file = Path(OUT_DIR, OUT_PREFIX + str(fp.stem) + '_' + solution_param).with_suffix('.out.txt')
    sys.stdout = open(out_file, 'w')
    nr3_vals = new_nr3_solution.get_data_frame(solution_param)
    dss_vals = dss_solution.get_data_frame(solution_param)
    test = ((nr3_vals - dss_vals).abs().max() <= tolerance).all()
    if test:
        print("TEST PASSED.")
    else:
        print("TEST FAILED.")
    compare_data_frames(nr3_vals, dss_vals, 'nr3', 'dss', solution_param)
    assert test
