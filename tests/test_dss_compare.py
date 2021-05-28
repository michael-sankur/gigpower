# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 23 March 2021
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
import csv

DSS_FILE_DIR = Path('./src/nr3_python/')
OUT_DIR_FBS = Path('./tests/test_compare_fbs_dss')
OUT_DIR_NR3 = Path('./tests/test_compare_nr3_dss')
REPORT_FILE = Path('./tests/test_report.csv')

# thresholds for passing tests
GENEROUS = .5
STRICT = .01


@pytest.mark.parametrize(
    "solution_param",
    [(param) for param in SolutionFBS.SOLUTION_PARAMS]
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
    "algorithm,out_dir",
    [
        ("NR3", OUT_DIR_NR3),
        ("FBS", OUT_DIR_FBS),
    ]
)
class TestDssCompare:

    summary_rows = [
        ['Algorithm', 'File', 'Zip', 'Parameter', 'Result', 'MaxErr']
    ]
    def setup_method(self):
        if not OUT_DIR_NR3.exists():
            os.makedirs(OUT_DIR_NR3)
        if not OUT_DIR_FBS.exists():
            os.makedirs(OUT_DIR_FBS)
    
    def teardown_method(cls):
        with open(REPORT_FILE, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(cls.summary_rows)

    def dss_solution(self, dss_file):
        fp = str(Path(DSS_FILE_DIR, dss_file))
        dss_solution = SolutionDSS(str(fp))
        dss_solution.solve()
        return dss_solution

    def new_fbs_solution(self, dss_file):
        fp = str(Path(DSS_FILE_DIR, dss_file))
        solution = SolutionFBS(fp)
        solution.solve()
        return solution

    def new_nr3_solution(self, dss_file):
        fp = str(Path(DSS_FILE_DIR, dss_file))
        solution = SolutionNR3(fp)
        solution.solve()
        return solution

    def test_dss_compare(self, algorithm, out_dir, zip_values,
                         zip_name, solution_param, dss_file, tolerance):
        """
        Compare the python FBS solution to the opendss solution.
        Writes output to OUTPUT FOLDER.
        """
        dss_solution = self.dss_solution(dss_file)
        if algorithm == 'NR3':
            new_solution = self.new_nr3_solution(dss_file)
        elif algorithm == 'FBS':
            new_solution = self.new_fbs_solution(dss_file)
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
            err = diffs.max()
            test = (err <= tolerance).all()
        else:
            if solution_param == 'Vmag':  # set Vmag tolerance to strict in all cases
                tolerance = STRICT
            err = (new_vals - dss_vals).abs().max()
            test = (err <= tolerance).all()
        if test:
            test_val = "Pass"
            print(f"TEST PASSED. {algorithm} CONVERGED = {new_solution.converged}")
        else:
            test_val = "Fail"
            print(f"TEST FAILED. {algorithm} CONVERGED = {new_solution.converged}")
        print(f"Zip values = {zip_values}")
        compare_data_frames(new_vals, dss_vals, algorithm, 'dss', solution_param)
        self.__class__.summary_rows.append([algorithm, dss_file, zip_name,
                                            solution_param, test_val, err])
        assert test
