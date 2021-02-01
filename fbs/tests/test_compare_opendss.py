# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: September 2020
# Compare solutions from opendss and powerflowpy.fbs

from fbs.utils import init_from_dss
from fbs.fbs import fbs, get_solution as get_fbs_solution
from fbs.dss_solve import solve_with_dss, getVMag

import pandas as pd
import pytest


# dss_files = ['fbs/tests/05n3ph_unbal/compare_opendss_05node_threephase_unbalanced_oscillation_03.dss']
# dss_files = ['fbs/tests/05n3ph_unbal/compare_opendss_05node_threephase_unbalanced_oscillation_03.dss',
#               'fbs/tests/06n3ph_unbal/06node_threephase_unbalanced.dss',
#               'fbs/tests/06n3ph_rad_unbal/06node_threephase_radial_unbalanced.dss']
# dss_files = ['fbs/tests/06n3ph_unbal/06node_threephase_unbalanced.dss']
# dss_files = ['fbs/tests/test_cases_dss/02node_threephase_unbalanced.dss']
# dss_files = ['/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/IEEE_13_Node/IEEE_13_Bus_allwye_noxfm_noreg.dss']
# dss_files = ['fbs/tests/IEEE_13_bus/IEEE_13_Bus_original.dss']
# dss_files = ['fbs/tests/IEEE_13_bus/IEEE_13_Bus_allwye.dss']
# dss_files = ['fbs/tests/IEEE_34_feeder_UB/IEEE_34_Bus_allwye.dss']
# dss_files = ['fbs/tests/IEEE_34_feeder_UB/IEEE_34_Bus_allwye_test01.dss']
# dss_files = ['fbs/tests/IEEE_34_feeder_UB/IEEE_34_Bus_allwye_test02.dss']
# dss_files = ['fbs/tests/IEEE_34_feeder_UB/IEEE_34_Bus_allwye_test03.dss']
# dss_files = ['fbs/tests/IEEE_34_feeder_UB/IEEE_34_Bus_allwye_test03a.dss']
# dss_files = ['/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/IEEE_34_feeder_UB/IEEE_34_Bus_allwye_noxfm_noreg.dss']
# dss_files = ['fbs/tests/IEEE_37_feeder_UB/IEEE_37_Bus_allwye.dss']
dss_files = ['/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/IEEE_37_feeder_UB/IEEE_37_Bus_allwye_noxfm_noreg.dss']


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
    dssV, dssI, dssStx, dssSrx, dssLoads, dss = dss_sol
    dssVMag = getVMag(dss)
    fbsVMag = fbs_sol.VMag_df()
    print(f"FBS iterations: {fbs_sol.iterations}\t FBS convergence:{fbs_sol.diff}\t FBS tolerance: {fbs_sol.tolerance}")
    V_maxDiff = compare_dfs(fbsV, dssV, "COMPARE V")
    VMag_maxDiff = compare_dfs(fbsVMag, dssVMag, "COMPARE VMag")
    I_maxDiff = compare_dfs(fbsI, dssI, "COMPARE I")
    Stx_maxDiff = compare_dfs(fbsStx, dssStx, "COMPARE Stx")
    Srx_maxDiff = compare_dfs(fbsSrx, dssSrx, "COMPARE Srx")

    fbsLoads = fbs_sV.multiply(1000)
    loads_maxDiff = compare_dfs(fbsLoads, dssLoads, "TOTAL LOAD POWERS")

    fbsLoads_sumPhase = fbsLoads.sum(axis = 0)
    dssLoads_sumPhase = dssLoads.sum(axis = 0)
    sumPhase_maxDiff = compare_dfs(fbsLoads_sumPhase, dssLoads_sumPhase,
    "TOTAL NODE POWERS - SUM OVER PHASE")

    assert (V_maxDiff <= tolerance).all()
    assert (VMag_maxDiff <= tolerance).all()
    assert (I_maxDiff <= tolerance).all()
    assert (Stx_maxDiff <= tolerance).all()
    assert (Srx_maxDiff <= tolerance).all()

    # convert node powers from kw
    assert (loads_maxDiff/1000 <= tolerance).all()


def compare_dfs(fbs_df: pd.DataFrame, dss_df: pd.DataFrame, title: str) -> None:
    """ helper method to compare fbs vs. dss and print comparisons """
    compare_cols = ['A.(fbs - dss)', 'B.(fbs - dss)', 'C.(fbs - dss)']
    compare = fbs_df.sub(dss_df)
    compare.columns = compare_cols
    pd.options.display.float_format = '{:.4f}'.format
    print(f"{title} - SUMMARY STATS")
    print(f"Max |fbs - dss|:\n{compare.abs().max()}\n")
    print(f"Sum |fbs - dss|:\n{compare.abs().sum()}\n")
    print(f"Avg |fbs - dss|:\n{compare.abs().mean()}\n")
    print(f"{title} - ERROR, |fbs - dss|:")
    print(compare, "\n")
    print(f"{title} - FBS RESULTS")
    print(fbs_df, "\n")
    print(f"{title} - DSS RESULTS")
    print(dss_df, "\n\n")
    return compare.abs().max()
