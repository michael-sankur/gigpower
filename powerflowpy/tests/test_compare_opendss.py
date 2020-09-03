# Compare the python FBS solution to the opendss solution.
import numpy as np
import opendssdirect as dss
from powerflowpy.utils import init_from_dss
from powerflowpy.fbs import *
import opendssdirect as dss

from math import tan, acos
import copy
import pandas as pd
import time
import re
import sys
import pytest

dss_file = 'powerflowpy/tests/05n3ph_unbal/compare_opendss_05node_threephase_unbalanced_oscillation_03.dss'

# construct the python FBS solution
def test_fbs_sol(dss_sol):
    #TODO: Figure out how to get Inode (iNR) and sV (sNR) at each node from dss and perform a compare
    fbs_sol= fbs(dss_file)
    network = init_from_dss(dss_file)
    fbsV, fbsI, fbsStx, fbsSrx = fbs_sol.V_df(), fbs_sol.I_df(), fbs_sol.Stx_df(), fbs_sol.Srx_df()
    dssV, dssI, dssStx, dssSrx = dss_sol
    print(f"FBS iterations: {fbs_sol.iterations}\t FBS convergence:{fbs_sol.diff}\t FBS tolerance: {fbs_sol.tolerance}")
    print(f"FBS loads: {network.to_dataframes()['Lines']}")
    print("\nCOMPARE V")
    compare_dfs(fbsV, dssV)
    print("\nCOMPARE I")
    compare_dfs(fbsI, dssI)
    print("\nCOMPARE Stx")
    compare_dfs(fbsStx, dssStx)
    print("\nCOMPARE Srx")
    compare_dfs(fbsSrx, dssSrx)
    print("Reprint fbsSrx")
    print(fbsSrx)

# construct the DSS solution. Copied form '20180601/opendss_nonvec_test_comparison.ipynb'

def compare_dfs(fbs_df : pd.DataFrame, dss_df : pd.DataFrame) -> None:
    """ helper method to compare fbs vs. dss and print comparisons """
    compare_cols = ['A.(fbs - dss)', 'B.(fbs - dss)', 'C.(fbs - dss)']
    compare = fbs_df.sub(dss_df)
    compare.columns = compare_cols
    concat = fbs_df.join(dss_df, lsuffix='.fbs', rsuffix='.dss')
    print("Max |diff|:")
    print(compare.abs().max())
    print(compare.join(concat))

@pytest.fixture
def dss_sol():
    """Run opendss's Solution.Solve on dss_file and save to a dictionary"""
    dss.run_command('Redirect ' + dss_file)
    # Set slack bus (sourcebus) voltage reference in p.u.
    SlackBusVoltage = 1.000
    dss.Vsources.PU(SlackBusVoltage)
    dss.Solution.Convergence(0.000000000001)

    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0
    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase

    # Solve power flow with OpenDSS file
    dss.Solution.Solve()
    if not dss.Solution.Converged:
        print('Initial Solution Not Converged. Check Model for Convergence')
    else:
        #Doing this solve command is required for GridPV, that is why the monitors
        #go under a reset process
        dss.Monitors.ResetAll()

        #set solution Params
        #setSolutionParams(dss,'daily',1,1,'off',1000000,30000)
        dss.Solution.Mode(1)
        dss.Solution.Number(1)
        dss.Solution.StepSize(1)
        dss.Solution.ControlMode(-1)
        dss.Solution.MaxControlIterations(1000000)
        dss.Solution.MaxIterations(30000)

        print(f'\nOpendss Iterations: {dss.Solution.Iterations()}\tOpendss Convergence: {dss.Solution.Convergence()}')

        VDSS = np.zeros((len(dss.Circuit.AllBusNames()), 3), dtype='complex')
        for k1 in range(len(dss.Circuit.AllBusNames())):
            dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k1])
            ph = np.asarray(dss.Bus.Nodes(), dtype='int')-1
            Vtemp = np.asarray(dss.Bus.PuVoltage())
            Vtemp = Vtemp[0:5:2] + 1j*Vtemp[1:6:2]
            VDSS[k1, ph] = Vtemp
        dssV = pd.DataFrame(VDSS, dss.Circuit.AllBusNames(), ['A', 'B', 'C'])

        IDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')
        for k1 in range(len(dss.Lines.AllNames())):
            dss.Lines.Name(dss.Lines.AllNames()[k1])

            ph = np.asarray(dss.CktElement.BusNames()[0].split('.')[1:], dtype='int')-1
            Imn = np.asarray(dss.CktElement.Currents())/Ibase
            Imn = Imn[0:int(len(Imn)/2)]
            Imn = Imn[0:5:2] + 1j*Imn[1:6:2]
            IDSS[k1,ph] = Imn
        dssI = pd.DataFrame(IDSS, dss.Lines.AllNames(), ['A', 'B', 'C'])

        STXDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')
        SRXDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')

        for k1 in range(len(dss.Lines.AllNames())):
            dss.Lines.Name(dss.Lines.AllNames()[k1])
            ph = np.asarray(dss.CktElement.BusNames()[0].split('.')[1:], dtype='int')-1
            Sk = np.asarray(dss.CktElement.Powers())/(Sbase/1000)
            STXtemp = Sk[0:int(len(Sk)/2)]
            SRXtemp = Sk[int(len(Sk)/2):]
            STXtemp = STXtemp[0:5:2] + 1j*STXtemp[1:6:2]
            SRXtemp = -(SRXtemp[0:5:2] + 1j*SRXtemp[1:6:2])
            STXDSS[k1, ph] = STXtemp
            SRXDSS[k1, ph] = SRXtemp

        dssStx = pd.DataFrame(STXDSS, dss.Lines.AllNames(), ['A', 'B', 'C'])
        dssSrx = pd.DataFrame(SRXDSS, dss.Lines.AllNames(), ['A', 'B', 'C'])
        print("OpenDSS Loads:")
        print(dss.CktElement.Powers())
        return dssV, dssI, dssStx, dssSrx
