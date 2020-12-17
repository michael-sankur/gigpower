# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: September 2020
# Compare solutions from opendss and powerflowpy.fbs

from fbs.utils import init_from_dss
from fbs.fbs import fbs, get_solution as get_fbs_solution
from fbs.dss_solve import solve_with_dss

import pandas as pd
import pytest


# dss_files = ['fbs/tests/05n3ph_unbal/compare_opendss_05node_threephase_unbalanced_oscillation_03.dss']
# dss_files = ['fbs/tests/05n3ph_unbal/compare_opendss_05node_threephase_unbalanced_oscillation_03.dss',
#               'fbs/tests/06n3ph_unbal/06node_threephase_unbalanced.dss',
#               'fbs/tests/06n3ph_rad_unbal/06node_threephase_radial_unbalanced.dss']
# dss_files = ['fbs/tests/06n3ph_unbal/06node_threephase_unbalanced.dss']
# dss_files = ['fbs/tests/test_cases_dss/02node_threephase_unbalanced.dss']
# dss_files = ['fbs/tests/IEEE_13_bus/IEEE_13_Bus_allwye_noxfm_noreg.dss']
# dss_files = ['fbs/tests/IEEE_13_bus/IEEE_13_Bus_original.dss']
# dss_files = ['fbs/tests/IEEE_13_bus/IEEE_13_Bus_allwye.dss']
# dss_files = ['fbs/tests/IEEE_34_feeder_UB/IEEE_34_Bus_allwye.dss']
# dss_files = ['fbs/tests/IEEE_37_feeder_UB/IEEE_37_Bus_allwye.dss']
dss_files = ['fbs/tests/IEEE_37_feeder_UB/IEEE_37_Bus_allwye_noxfm_noreg.dss']



def test_all():
    """
    Compare the python FBS solution to the opendss solution.
    To save output, run from command line:
    pytest ./fbs/tests/test_compare_opendss.py -s >[OUTPUT FILE PATH]
    Examples:
    $ pytest ./fbs/tests/test_compare_opendss.py -s > ./fbs/tests/06n3ph_unbal/06n3ph_out.txt
    $ pytest ./fbs/tests/test_compare_opendss.py -s > ./fbs/tests/05n3ph_unbal/05n3ph_out.txt
    """
    tolerance = 10**-1
    for file in dss_files:
        print(f'\n\nTesting File: {file}')
        compare_fbs_sol(file, tolerance)


# construct the python FBS solution
def compare_fbs_sol(dss_file, tolerance):
    # TODO: Figure out how to get Inode (iNR) and sV (sNR) at each node from dss and perform a compare
    # solve with fbs
    network = init_from_dss(dss_file)
    fbs_sol = get_fbs_solution(fbs(network))

    # solve with dss
    dss_sol = solve_with_dss(dss_file)

    fbsV, fbsI, fbsStx, fbsSrx, fbs_sV = fbs_sol.V_df(), fbs_sol.I_df(), fbs_sol.Stx_df(), \
        fbs_sol.Srx_df(), fbs_sol.sV_df()
    dssV, dssI, dssStx, dssSrx, dssLoads = dss_sol
    print(f"FBS iterations: {fbs_sol.iterations}\t FBS convergence:{fbs_sol.diff}\t FBS tolerance: {fbs_sol.tolerance}")
    print("\n\nCOMPARE V")
    V_maxDiff = compare_dfs(fbsV, dssV)
    print("\n\nCOMPARE I")
    I_maxDiff = compare_dfs(fbsI, dssI)
    print("\n\nCOMPARE Stx")
    Stx_maxDiff = compare_dfs(fbsStx, dssStx)
    print("\n\nCOMPARE Srx")
    Srx_maxDiff = compare_dfs(fbsSrx, dssSrx)

    print("\n\nCOMPARE TOTAL NODE POWERS\n")
    fbsLoads = fbs_sV.multiply(1000)

    loads_maxDiff = compare_dfs(fbsLoads, dssLoads)

    assert (V_maxDiff <= tolerance).all()
    assert (I_maxDiff <= tolerance).all()
    assert (Stx_maxDiff <= tolerance).all()
    assert (Srx_maxDiff <= tolerance).all()

    # convert node powers from kw
    assert (loads_maxDiff/1000 <= tolerance).all()


# construct the DSS solution. Copied form '20180601/opendss_nonvec_test_comparison.ipynb'
def compare_dfs(fbs_df: pd.DataFrame, dss_df: pd.DataFrame) -> None:
    """ helper method to compare fbs vs. dss and print comparisons """
    compare_cols = ['A.(fbs - dss)', 'B.(fbs - dss)', 'C.(fbs - dss)']
    compare = fbs_df.sub(dss_df)
    compare.columns = compare_cols
    pd.options.display.float_format = '{:.4f}'.format
    print("Max |diff|:")
    print(compare.abs().max())
    print(compare, "\n")
    print("FBS results")
    print(fbs_df, "\n")
    print("DSS results")
    print(dss_df)
    return compare.abs().max()
