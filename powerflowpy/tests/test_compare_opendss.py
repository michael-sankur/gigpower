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
def test_fbs_sol():
    solution = fbs(dss_file)
    print("\nPython FBS Solution")
    solution.print_solution()
    print()
# construct the DSS solution. Copied form '20180601/opendss_nonvec_test_comparison.ipynb'
# @pytest.fixture
def test_dss_sol():
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
        #Easy process to get all names and count of loads, a trick to avoid
        #some more lines of code
        TotalLoads = dss.Loads.Count()
        AllLoadNames = dss.Loads.AllNames()
        print('OpenDSS Model Compliation Done.')
        print('Iterations: ', dss.Solution.Iterations())
        print('Tolerance: ', dss.Solution.Convergence())

        print('')

        # print(dss.Solution.Converged())

        # Print number of buses, and bus names
        print(len(dss.Circuit.AllBusNames()))
        print(dss.Circuit.AllBusNames())

        # Print number of loads, and load names
        print(len(dss.Loads.AllNames()))
        print(dss.Loads.AllNames())

        print('')

        VDSS = np.zeros((3, len(dss.Circuit.AllBusNames())), dtype='complex')

        for k1 in range(len(dss.Circuit.AllBusNames())):
            dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k1])
        #     print(dss.Circuit.AllBusNames()[k1])
        #     print(dss.Bus.Nodes())

        #     print('puVOTLAGES - LN CARTESIAN')
        #     print(dss.Bus.PuVoltage())

            ph = np.asarray(dss.Bus.Nodes(), dtype='int')-1

            Vtemp = np.asarray(dss.Bus.PuVoltage())
            Vtemp = Vtemp[0:5:2] + 1j*Vtemp[1:6:2]

        #     print(np.asarray(dss.Bus.Nodes(),'int'))

            VDSS[ph, k1] = Vtemp

        #     VDSS[np.asarray(dss.Bus.Nodes(),'int'),k1] = np.array(dss.Bus.PuVoltage()[0:5:2] + 1j*dss.Bus.PuVoltage()[1:6:2])


        #     VDSS[dss.Bus.Nodes()-1,k1] = dss.Bus.PuVoltage()[0:2:5]
        #     for k2 in range(len(dss.Bus.Nodes())):
        #         VDSS[int(dss.Bus.Nodes()[k2])-1,k1] = dss.Bus.PuVoltage()[2*k2] + 1j*dss.Bus.PuVoltage()[2*k2+1]

        print('VDSS\n', np.round(VDSS, decimals=6))


        IDSS = np.zeros((3, len(dss.Lines.AllNames())), dtype='complex')

        for k1 in range(len(dss.Lines.AllNames())):
            dss.Lines.Name(dss.Lines.AllNames()[k1])
        #     print(dss.Lines.AllNames()[k1])
            ph = np.asarray(dss.CktElement.BusNames()[0].split('.')[1:], dtype='int')-1
            Imn = np.asarray(dss.CktElement.Currents())/Ibase
        #     print(Imn)
            Imn = Imn[0:int(len(Imn)/2)]
        #     print(Imn)
            Imn = Imn[0:5:2] + 1j*Imn[1:6:2]
        #     print(Imn)
            IDSS[ph, k1] = Imn
        #     print('')

        print('IDSS\n', np.round(IDSS, decimals=6))

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
        assert True
