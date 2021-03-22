# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: February 2021
# Compare node powers from opendss and LinDistFlow.fbs in order to
# to highlight discrepancies with dss.CktElement.Powers()

from fbs.utils import init_from_dss
from fbs.fbs import fbs, get_solution as get_fbs_solution
from fbs import dss_solve
from fbs.dss_solve import solve_with_dss, getTotalNodePowers, getNominalNodePower, \
    getTotalNodePowers_CktElement

import pandas as pd

CUR_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
dss_file = CUR_DIR + 'IEEE_13_Node/IEEE_13_Bus_allwye_noxfm_noreg.dss'
# dss_file = 'IEEE_34_feeder_UB/IEEE_34_Bus_allwye_noxfm_noreg.dss'
# dss_file = 'IEEE_37_Bus_allwye_noxfm_noreg.dss'


def relabelColumns(df, prefix):
    cols = df.columns
    return df.rename(columns={ph: f"{prefix}_{ph}" for ph in cols})


def compare_nominals(fbs_sol, dss):
    """
    param fbs_sol: A solved FBS Solution object
    param dss: a solved DSS object
    """
    dssNominal = getNominalNodePower(dss)
    fbsNominal = fbs_sol.nomNodePwrs_df() * 1000
    dfs = [dssNominal, fbsNominal]
    prefixes = ['dssNom', 'fbsNom']
    renamed = [relabelColumns(df, prefix) for df, prefix in zip(dfs, prefixes)]
    result = pd.concat(renamed, axis=1)
    result = result[['dssNom_A', 'fbsNom_A',
                     'dssNom_B', 'fbsNom_B',
                     'dssNom_C', 'fbsNom_C']]
    print("\nNOMINAL POWERS, dss, fbs")
    print(result, '\n')
    print("NOMINAL POWERS, dss - fbs")
    print(dssNominal - fbsNominal)
    print('~' * 100, '\n')
    return dssNominal - fbsNominal


def test_node_power_compare():
    """
    Compare the python FBS solution to the opendss solution.
    To save output, run from command line:
    pytest ./fbs/tests/test_node_power_compare.py -s > [OUTPUT FILE PATH]
    Examples:
    $ pytest ./fbs/tests/test_node_power_compare.py -s > fbs/tests/IEEE_13_bus/IEEE_13_Bus_allwye_fbs_out.txt
    $
    """
    pd.options.display.float_format = '{:.1f}'.format
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 150)
    print(f'\n\nTesting File: {dss_file}')

    network = init_from_dss(dss_file)
    fbs_sol = get_fbs_solution(fbs(network))
    dssV, dssI, dssStx, dssSrx, dssNodePowers, dss = solve_with_dss(dss_file)

    # NOMINALS ------------------------------------------------------
    compare_nominals(fbs_sol, dss)

    # LOAD POWERS BY BUS --------------------------------------------
    fbsLoadPowers = fbs_sol.getLoadPowers() * 1000

    dssCktElementLoads = dss_solve.getLoadPowers_CktElement(dss)
    dssCalcLoads = dss_solve.getLoadPowers(dss) * 1000
    dfs = [dssCktElementLoads, dssCalcLoads, fbsLoadPowers]
    prefixes = ['dss1', 'dss2', 'fbs']
    renamed = [relabelColumns(df, prefix) for df, prefix in zip(dfs, prefixes)]
    result = pd.concat(renamed, axis=1)
    result = result[['dss1_A', 'dss2_A', 'fbs_A',
                     'dss1_B', 'dss2_B', 'fbs_B',
                     'dss1_C', 'dss2_C', 'fbs_C']]

    print('LOAD POWERS BY BUS: dss1, dss2, fbs')
    print("dss1: Values from summing over dss.CktElement.Powers()")
    print("dss2: Calculated from opendss voltage solution and zip values\n")
    print(result, '\n')

    dss_diff = relabelColumns(dssCktElementLoads - dssCalcLoads, 'dss1 - dss2')
    dssCktElement_v_fbs = relabelColumns(dssCktElementLoads - fbsLoadPowers, 'dss1 - fbs')
    dssCalc_v_fbs = relabelColumns(dssCalcLoads - fbsLoadPowers, 'dss2 - fbs')
    print('LOAD POWERS BY BUS: dss1 - dss2')
    print(dss_diff, '\n')
    print('LOAD POWERS BY BUS: dss1 - fbs')
    print(dssCktElement_v_fbs, '\n')
    print('LOAD POWERS BY BUS: dss2 - fbs')
    print(dssCalc_v_fbs, '\n')

    # CAP POWERS BY BUS --------------------------------------------
    fbsCapPowers = fbs_sol.getCapPowers() * 1000
    dssCktElementCaps = dss_solve.getCapacitorPowers_CktElement(dss)
    dssCalcCaps = dss_solve.getCapacitorPowers(dss) * 1000
    dfs = [dssCktElementCaps, dssCalcCaps, fbsCapPowers]
    prefixes = ['dss1', 'dss2', 'fbs']
    renamed = [relabelColumns(df, prefix) for df, prefix in zip(dfs, prefixes)]
    result = pd.concat(renamed, axis=1)
    result = result[['dss1_A', 'dss2_A', 'fbs_A',
                     'dss1_B', 'dss2_B', 'fbs_B',
                     'dss1_C', 'dss2_C', 'fbs_C']]

    print('CAP POWERS BY BUS: dss1, dss2, fbs')
    print("dss1: Values from summing over dss.CktElement.Powers()")
    print("dss2: Calculated from opendss voltage solution and zip values\n")

    print(result, '\n')

    dss_diff = relabelColumns(dssCktElementCaps - dssCalcCaps, 'dss1 - dss2')
    dssCktElement_v_fbs = relabelColumns(dssCktElementCaps - fbsCapPowers, 'dss1 - fbs')
    dssCalc_v_fbs = relabelColumns(dssCalcCaps - fbsCapPowers, 'dss2 - fbs')
    print('CAP POWERS BY BUS: dss1 - dss2')
    print(dss_diff, '\n')
    print('CAP POWERS BY BUS: dss1 - fbs')
    print(dssCktElement_v_fbs, '\n')
    print('CAP POWERS BY BUS: dss2 - fbs')
    print(dssCalc_v_fbs, '\n')

    # TOTAL POWERS BY BUS------------------------------------------------
    fbsNodePowers = fbs_sol.sV_df() * 1000
    dssCktElementPowers = getTotalNodePowers_CktElement(dss)
    dssCalcNodePowers = getTotalNodePowers(dss) * 1000
    fbsNodePowers = fbs_sol.sV_df() * 1000
    dfs = [dssCktElementPowers, dssCalcNodePowers, fbsNodePowers]
    prefixes = ['dss1', 'dss2', 'fbs']
    renamed = [relabelColumns(df, prefix) for df, prefix in zip(dfs, prefixes)]
    result = pd.concat(renamed, axis=1)
    result = result[['dss1_A', 'dss2_A', 'fbs_A',
                     'dss1_B', 'dss2_B', 'fbs_B',
                     'dss1_C', 'dss2_C', 'fbs_C']]

    print('TOTAL BUS POWERS: dss1, dss2, fbs')
    print("dss1: Values from summing over dss.CktElement.Powers()")
    print("dss2: Calculated from opendss voltage solution and zip values\n")

    print(result, '\n')

    dss_diff = relabelColumns(dssCktElementPowers - dssCalcNodePowers, 'dss1 - dss2')
    dssCktElement_v_fbs = relabelColumns(dssCktElementPowers - fbsNodePowers, 'dss1 - fbs')
    dssCalc_v_fbs = relabelColumns(dssCalcNodePowers - fbsNodePowers, 'dss2 - fbs')
    print('TOTAL BUS POWERS: dss1 - dss2')
    print(dss_diff, '\n')
    print('TOTAL BUS POWERS: dss1 - fbs')
    print(dssCktElement_v_fbs, '\n')
    print('TOTAL BUS POWERS: dss2 - fbs')
    print(dssCalc_v_fbs, '\n')
