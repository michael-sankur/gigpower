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
    fbs_sol= fbs(dss_file)
    network = init_from_dss(dss_file)
    fbsV, fbsI = fbs_sol.V_df(), fbs_sol.I_df()
    dssV, dssI = dss_sol
    compare_cols = ['A.(fbs - dss)', 'B.(fbs - dss)', 'C.(fbs - dss)']
    print("\nCOMPARE V")
    compareV = fbsV.sub(dssV)
    compareV.columns = compare_cols
    concatV = dssV.join(fbsV, lsuffix='.dss', rsuffix='.fbs')
    print("Max |diff|:")
    print(compareV.abs().max())
    print(compareV.join(concatV))

    print("\nCOMPARE I")
    compareI = fbsI.sub(dssI)
    compareI.columns = compare_cols
    concatI = dssI.join(fbsI, lsuffix='.dss', rsuffix='.fbs')
    print("Max |diff|:")
    print(compareI.abs().max())
    print(compareI.join(concatI))
# construct the DSS solution. Copied form '20180601/opendss_nonvec_test_comparison.ipynb'

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
        print('Initial Model Converged. Proceeding to Next Step.')
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

        print('Opendss Iterations: ', dss.Solution.Iterations())
        print('Opendss Tolerance: ', dss.Solution.Convergence())

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

        STXDSS = np.zeros((3, len(dss.Lines.AllNames())), dtype='complex')
        SRXDSS = np.zeros((3, len(dss.Lines.AllNames())), dtype='complex')

        for k1 in range(len(dss.Lines.AllNames())):
            dss.Lines.Name(dss.Lines.AllNames()[k1])
        #     print(dss.Lines.AllNames()[k1])
        #     print(dss.CktElement.BusNames())
            ph = np.asarray(dss.CktElement.BusNames()[0].split('.')[1:], dtype='int')-1
        #     print(ph)
            Sk = np.asarray(dss.CktElement.Powers())/(Sbase/1000)
        #     print(Sk)

        #     print(Sk[0:int(len(Sk)/2)])
        #     print(Sk[int(len(Sk)/2):])

            STXtemp = Sk[0:int(len(Sk)/2)]
            SRXtemp = Sk[int(len(Sk)/2):]

            STXtemp = STXtemp[0:5:2] + 1j*STXtemp[1:6:2]
        #     print(STXtemp)

            SRXtemp = -(SRXtemp[0:5:2] + 1j*SRXtemp[1:6:2])
        #     print(SRXtemp)

            STXDSS[ph, k1] = STXtemp
            SRXDSS[ph, k1] = SRXtemp

        print('STXDSS\n', np.round(STXDSS, decimals=6))
        print('SRXDSS\n', np.round(SRXDSS, decimals=6))

        print('|VDSS|\n', np.round(np.abs(VDSS), decimals=6))
        print('<VDSS\n', np.round(180/np.pi*np.angle(VDSS), decimals=6))
        print('D<VDSS\n', 180/np.pi*np.angle(VDSS) - 180/np.pi*np.angle(VDSS[:, [0]]))
        return dssV, dssI
