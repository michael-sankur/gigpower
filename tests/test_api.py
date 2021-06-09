# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 28 May 2021
# Smell tests to verify Solution API functions

from circuit_mapper.solution import Solution
from circuit_mapper.solution_dss import SolutionDSS
from circuit_mapper.solution_fbs import SolutionFBS
from circuit_mapper.solution_nr3 import SolutionNR3
from circuit_mapper.pretty_print import compare_data_frames

from circuit_mapper.utils import get_nominal_bus_powers

import pytest
from pathlib import Path
import opendssdirect as dss
import pandas as pd

DSS_FILE_DIR = Path('./src/nr3_python/')

@pytest.mark.parametrize(
    "dss_file",
    [
        ('IEEE_13_Bus_allwye.dss'),
        ('IEEE_13_Bus_allwye_noxfm_noreg.dss'),
        ('IEEE_34_Bus_allwye.dss'),
        ('IEEE_34_Bus_allwye_noxfm_noreg.dss'),
        ('IEEE_37_Bus_allwye.dss'),
        ('IEEE_37_Bus_allwye_noxfm_noreg.dss')
    ]
)
@pytest.mark.parametrize(
    "algorithm",
    [
        (SolutionNR3),
        (SolutionFBS),
        (SolutionDSS)
    ]
)
class TestSolutionDFs:

    def get_solution(self, dss_file, algorithm):
        fp = str(Path(DSS_FILE_DIR, dss_file))
        solution = algorithm(str(fp))
        solution.solve()
        return solution

    def test_dfs(self, dss_file, algorithm):
        """
        Run calls to get Solution.V, Solution.I, Solution.sV, Solution.VMag
        as data frames
        """
        solution = self.get_solution(dss_file, algorithm)
        for param in Solution.SOLUTION_PARAMS:
            df = solution.get_data_frame(param)
            pytest.assume(not(df.empty))  # make sure df is not empty

    def test_dfs_orient(self, dss_file, algorithm):
        """
        Run calls to get Solution.V, Solution.I, Solution.sV, Solution.VMag
        as data frames with both orientations (rows, columns) and make sure
        that they have transposed shapes
        """
        solution = self.get_solution(dss_file, algorithm)
        for param in Solution.SOLUTION_PARAMS:
            df_rows = solution.get_data_frame(param, orient='rows')
            df_cols = solution.get_data_frame(param, orient='cols')
            pytest.assume(df_rows.shape[-1::-1] == df_cols.shape)
            # check that 3 phases are oriented correctly
            pytest.assume(df_rows.shape[1] == 3)
            pytest.assume(df_cols.shape[0] == 3)

    def test_nominals(self, dss_file, algorithm):
        """
        Make sure that Circuit class's nominal powers match those from 
        opendss' api
        """
        solution = self.get_solution(dss_file, algorithm)
        solution_nominals = solution.get_nominal_bus_powers(orient='rows')

        # get a fresh dss object for each new dss file
        fp = str(Path(DSS_FILE_DIR, dss_file))
        dss.run_command('Redirect ' + fp)
        dss.Solution.Solve()
        dss_nominals = get_nominal_bus_powers(dss)

        pd.testing.assert_frame_equal(solution_nominals, dss_nominals)