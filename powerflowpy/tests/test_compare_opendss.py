# Compare the python FBS solution to the opendss solution.
# To save output, run from command line:
# pytest ./powerflowpy/tests/test_compare_opendss.py -s >[OUTPUT FILE PATH]
# Ex:
# pytest ./powerflowpy/tests/test_compare_opendss.py -s > ./powerflowpy/tests/06n3ph_unbal/06n3ph_out.txt
# pytest ./powerflowpy/tests/test_compare_opendss.py -s > ./powerflowpy/tests/05n3ph_unbal/05n3ph_out.txt

import numpy as np
import opendssdirect as dss
from powerflowpy.utils import init_from_dss
from powerflowpy.fbs import *
from powerflowpy.dss_solve import solve_with_dss
import opendssdirect as dss

from math import tan, acos
import copy
import pandas as pd
import time
import re
import sys
import pytest

# dss_files = ['powerflowpy/tests/05n3ph_unbal/compare_opendss_05node_threephase_unbalanced_oscillation_03.dss', 'powerflowpy/tests/06n3ph_unbal/06node_threephase_unbalanced.dss', 'powerflowpy/tests/06n3ph_rad_unbal/06node_threephase_radial_unbalanced.dss']
dss_files = ['powerflowpy/tests/06n3ph_rad_unbal/06node_threephase_radial_unbalanced.dss']

def test_all():
    tolerance = .01
    for file in dss_files:
        print(f'\n\nTesting File: {file}')
        compare_fbs_sol(file, tolerance)

# construct the python FBS solution
def compare_fbs_sol(dss_file, tolerance):
    #TODO: Figure out how to get Inode (iNR) and sV (sNR) at each node from dss and perform a compare
    network = init_from_dss(dss_file)
    fbs_sol = fbs(network, False)
    dss_sol = solve_with_dss(dss_file)

    fbsV, fbsI, fbsStx, fbsSrx = fbs_sol.V_df(), fbs_sol.I_df(), fbs_sol.Stx_df(), fbs_sol.Srx_df()
    dssV, dssI, dssStx, dssSrx = dss_sol
    print(f"FBS iterations: {fbs_sol.iterations}\t FBS convergence:{fbs_sol.diff}\t FBS tolerance: {fbs_sol.tolerance}")
    print("\nCOMPARE V")
    V_maxDiff = compare_dfs(fbsV, dssV)
    print("\nCOMPARE I")
    I_maxDiff = compare_dfs(fbsI, dssI)
    print("\nCOMPARE Stx")
    Stx_maxDiff = compare_dfs(fbsStx, dssStx)
    print("\nCOMPARE Srx")
    Srx_maxDiff = compare_dfs(fbsSrx, dssSrx)

    assert (V_maxDiff <= tolerance).all()
    assert (I_maxDiff <= tolerance).all()
    assert (Stx_maxDiff <= tolerance).all()
    assert (Srx_maxDiff <= tolerance).all()

# construct the DSS solution. Copied form '20180601/opendss_nonvec_test_comparison.ipynb'

def compare_dfs(fbs_df : pd.DataFrame, dss_df : pd.DataFrame) -> None:
    """ helper method to compare fbs vs. dss and print comparisons """
    compare_cols = ['A.(fbs - dss)', 'B.(fbs - dss)', 'C.(fbs - dss)']
    compare = fbs_df.sub(dss_df)
    compare.columns = compare_cols
    concat = fbs_df.join(dss_df, lsuffix='.fbs', rsuffix='.dss')
    print("Max |diff|:")
    print(compare.abs().max())
    print(compare)
    print("FBS results")
    print(fbs_df)
    print("DSS results")
    print(dss_df)
    return compare.abs().max()
