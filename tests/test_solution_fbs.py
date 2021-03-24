# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 21 March 2021
# Test refactor of nr3 module

import numpy as np
import pytest
import pandas as pd

from circuit_mapper.circuit import Circuit
from circuit_mapper.solution_fbs import SolutionFBS

from fbs.fbs.utils import init_from_dss
from fbs.fbs.fbs import fbs, get_solution as get_fbs_solution

LOCAL_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
DSS_FILE = LOCAL_DIR + 'IEEE_13_Node/IEEE_13_Bus_allwye_noxfm_noreg.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])

TOLERANCE = 10e-2

@pytest.fixture
def circuit():
    """ map Circuit object once for use in all tests """
    return Circuit(DSS_FILE)


@pytest.fixture
def old_fbs_solution():
    network = init_from_dss(DSS_FILE)
    return get_fbs_solution(fbs(network, 100))


@pytest.fixture
def new_fbs_solution():
    fbs = SolutionFBS(DSS_FILE)
    fbs.maxiter = 100
    fbs.solve()
    return fbs


def test_fbs_V(old_fbs_solution, new_fbs_solution):
    new = new_fbs_solution.get_data_frame('V')
    old = old_fbs_solution.V_df()
    print(new)
    print(old)
    assert ((new - old).abs().max() <= TOLERANCE).all()
    # assert new.equals(old)
    
def test_fbs_I(old_fbs_solution, new_fbs_solution):
    new = new_fbs_solution.get_data_frame('I')
    old = old_fbs_solution.I_df()
    print(new)
    print(old)
    assert ((new - old).abs().max() <= TOLERANCE).all()
    # assert new.equals(old)

# def test_fbs_Stx(old_fbs_solution, new_fbs_solution):
#     assert (old_fbs_solution.Stx_df() == new_fbs_solution.Stx_df()).all()


# def test_fbs_Srx(old_fbs_solution, new_fbs_solution):
#     assert (old_fbs_solution.Srx_df() == new_fbs_solution.Srx_df()).all()


def test_fbs_sV(old_fbs_solution, new_fbs_solution):
    new = new_fbs_solution.get_data_frame('sV')
    old = old_fbs_solution.sV_df()
    print(new)
    print(old)
    assert ((new - old).abs().max() <= TOLERANCE).all()
    # assert new.equals(old)
    
# def test_powers_nominal(old_fbs_solution, new_fbs_solution):
#     old_fbsNominal = old_fbs_solution.nomNodePwrs_df() * 1000
#     new_fbsNominal = new_fbs_solution.get_nominal_powers()


# def test_node_power_compare():
#     """
#     Compare the python FBS solution to the opendss solution.
#     To save output, run from command line:
#     pytest ./fbs/tests/test_node_power_compare.py -s > [OUTPUT FILE PATH]
#     Examples:
#     $ pytest ./fbs/tests/test_node_power_compare.py -s > fbs/tests/IEEE_13_bus/IEEE_13_Bus_allwye_fbs_out.txt
#     $
#     """
#     pd.options.display.float_format = '{:.1f}'.format
#     pd.set_option('display.max_rows', 500)
#     pd.set_option('display.max_columns', 500)
#     pd.set_option('display.width', 150)
#     print(f'\n\nTesting File: {dss_file}')

#     network = init_from_dss(dss_file)
#     fbs_sol = get_fbs_solution(fbs(network))
#     dssV, dssI, dssStx, dssSrx, dssNodePowers, dss = solve_with_dss(dss_file)

#     # NOMINALS ------------------------------------------------------
#     compare_nominals(fbs_sol, dss)

#     # LOAD POWERS BY BUS --------------------------------------------
#     fbsLoadPowers = fbs_sol.getLoadPowers() * 1000

#     dssCktElementLoads = dss_solve.getLoadPowers_CktElement(dss)
#     dssCalcLoads = dss_solve.getLoadPowers(dss) * 1000
#     dfs = [dssCktElementLoads, dssCalcLoads, fbsLoadPowers]
#     prefixes = ['dss1', 'dss2', 'fbs']
#     renamed = [relabelColumns(df, prefix) for df, prefix in zip(dfs, prefixes)]
#     result = pd.concat(renamed, axis=1)
#     result = result[['dss1_A', 'dss2_A', 'fbs_A',
#                      'dss1_B', 'dss2_B', 'fbs_B',
#                      'dss1_C', 'dss2_C', 'fbs_C']]

#     print('LOAD POWERS BY BUS: dss1, dss2, fbs')
#     print("dss1: Values from summing over dss.CktElement.Powers()")
#     print("dss2: Calculated from opendss voltage solution and zip values\n")
#     print(result, '\n')

#     dss_diff = relabelColumns(dssCktElementLoads - dssCalcLoads, 'dss1 - dss2')
#     dssCktElement_v_fbs = relabelColumns(dssCktElementLoads - fbsLoadPowers, 'dss1 - fbs')
#     dssCalc_v_fbs = relabelColumns(dssCalcLoads - fbsLoadPowers, 'dss2 - fbs')
#     print('LOAD POWERS BY BUS: dss1 - dss2')
#     print(dss_diff, '\n')
#     print('LOAD POWERS BY BUS: dss1 - fbs')
#     print(dssCktElement_v_fbs, '\n')
#     print('LOAD POWERS BY BUS: dss2 - fbs')
#     print(dssCalc_v_fbs, '\n')

#     # CAP POWERS BY BUS --------------------------------------------
#     fbsCapPowers = fbs_sol.getCapPowers() * 1000
#     dssCktElementCaps = dss_solve.getCapacitorPowers_CktElement(dss)
#     dssCalcCaps = dss_solve.getCapacitorPowers(dss) * 1000
#     dfs = [dssCktElementCaps, dssCalcCaps, fbsCapPowers]
#     prefixes = ['dss1', 'dss2', 'fbs']
#     renamed = [relabelColumns(df, prefix) for df, prefix in zip(dfs, prefixes)]
#     result = pd.concat(renamed, axis=1)
#     result = result[['dss1_A', 'dss2_A', 'fbs_A',
#                      'dss1_B', 'dss2_B', 'fbs_B',
#                      'dss1_C', 'dss2_C', 'fbs_C']]

#     print('CAP POWERS BY BUS: dss1, dss2, fbs')
#     print("dss1: Values from summing over dss.CktElement.Powers()")
#     print("dss2: Calculated from opendss voltage solution and zip values\n")

#     print(result, '\n')

#     dss_diff = relabelColumns(dssCktElementCaps - dssCalcCaps, 'dss1 - dss2')
#     dssCktElement_v_fbs = relabelColumns(dssCktElementCaps - fbsCapPowers, 'dss1 - fbs')
#     dssCalc_v_fbs = relabelColumns(dssCalcCaps - fbsCapPowers, 'dss2 - fbs')
#     print('CAP POWERS BY BUS: dss1 - dss2')
#     print(dss_diff, '\n')
#     print('CAP POWERS BY BUS: dss1 - fbs')
#     print(dssCktElement_v_fbs, '\n')
#     print('CAP POWERS BY BUS: dss2 - fbs')
#     print(dssCalc_v_fbs, '\n')

#     # TOTAL POWERS BY BUS------------------------------------------------
#     fbsNodePowers = fbs_sol.sV_df() * 1000
#     dssCktElementPowers = getTotalNodePowers_CktElement(dss)
#     dssCalcNodePowers = getTotalNodePowers(dss) * 1000
#     fbsNodePowers = fbs_sol.sV_df() * 1000
#     dfs = [dssCktElementPowers, dssCalcNodePowers, fbsNodePowers]
#     prefixes = ['dss1', 'dss2', 'fbs']
#     renamed = [relabelColumns(df, prefix) for df, prefix in zip(dfs, prefixes)]
#     result = pd.concat(renamed, axis=1)
#     result = result[['dss1_A', 'dss2_A', 'fbs_A',
#                      'dss1_B', 'dss2_B', 'fbs_B',
#                      'dss1_C', 'dss2_C', 'fbs_C']]

#     print('TOTAL BUS POWERS: dss1, dss2, fbs')
#     print("dss1: Values from summing over dss.CktElement.Powers()")
#     print("dss2: Calculated from opendss voltage solution and zip values\n")

#     print(result, '\n')

#     dss_diff = relabelColumns(dssCktElementPowers - dssCalcNodePowers, 'dss1 - dss2')
#     dssCktElement_v_fbs = relabelColumns(dssCktElementPowers - fbsNodePowers, 'dss1 - fbs')
#     dssCalc_v_fbs = relabelColumns(dssCalcNodePowers - fbsNodePowers, 'dss2 - fbs')
#     print('TOTAL BUS POWERS: dss1 - dss2')
#     print(dss_diff, '\n')
#     print('TOTAL BUS POWERS: dss1 - fbs')
#     print(dssCktElement_v_fbs, '\n')
#     print('TOTAL BUS POWERS: dss2 - fbs')
#     print(dssCalc_v_fbs, '\n')
