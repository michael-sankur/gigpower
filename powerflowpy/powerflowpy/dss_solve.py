# # Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 10 September 2020
# run opendss' solver on a dss file. Copied form '20180601/opendss_nonvec_test_comparison.ipynb'

import opendssdirect as dss
import numpy as np
import pandas as pd

def solve_with_dss(dss_file):
    dss.run_command('Redirect ' + dss_file)

    # Set slack bus (sourcebus) voltage reference in p.u.
    SlackBusVoltage = 1.000
    dss.Vsources.PU(SlackBusVoltage)
    dss.Solution.Convergence(0.000000000001)

    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0
    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase


    # Solve power flow with OpenDSS file, and updating loads at each timestep.
    # Code based on https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Version8/Distrib/Examples/Python/Python-to-OpenDSS%20Control%20Interface.pdf
    dss.Solution.MaxControlIterations(1000000)
    originalSteps = dss.Solution.Number()
    dss.Solution.Number(1)
    # save original loads
    orig_loads_data = dss.utils.loads_to_dataframe()
    orig_loads_data = orig_loads_data.transpose()

    for stepNumber in range(originalSteps): # simulate timesteps

        # set loads
        for load_name in dss.Loads.AllNames():
            load_data = orig_loads_data[load_name]
            dss.Loads.Name(load_name)
            dss.Loads.kW( load_data['kW'] )
            dss.Loads.kvar( load_data['kvar'] )

        #run solve for this timestep
        dss.Solution.SolveSnap()
        dss.Solution.FinishTimeStep()

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
        # dssV = pd.DataFrame(VDSS, dss.Circuit.AllBusNames(), ['A', 'B', 'C'])

        # IDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')
        # for k1 in range(len(dss.Lines.AllNames())):
        #     dss.Lines.Name(dss.Lines.AllNames()[k1])

        #     ph = np.asarray(dss.CktElement.BusNames()[0].split('.')[1:], dtype='int')-1
        #     Imn = np.asarray(dss.CktElement.Currents())/Ibase
        #     Imn = Imn[0:int(len(Imn)/2)]
        #     Imn = Imn[0:5:2] + 1j*Imn[1:6:2]
        #     IDSS[k1,ph] = Imn
        # dssI = pd.DataFrame(IDSS, dss.Lines.AllNames(), ['A', 'B', 'C'])

        # STXDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')
        # SRXDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')

        # for k1 in range(len(dss.Lines.AllNames())):
        #     dss.Lines.Name(dss.Lines.AllNames()[k1])
        #     ph = np.asarray(dss.CktElement.BusNames()[0].split('.')[1:], dtype='int')-1
        #     Sk = np.asarray(dss.CktElement.Powers())/(Sbase/1000)
        #     STXtemp = Sk[0:int(len(Sk)/2)]
        #     SRXtemp = Sk[int(len(Sk)/2):]
        #     STXtemp = STXtemp[0:5:2] + 1j*STXtemp[1:6:2]
        #     SRXtemp = -(SRXtemp[0:5:2] + 1j*SRXtemp[1:6:2])
        #     STXDSS[k1, ph] = STXtemp
        #     SRXDSS[k1, ph] = SRXtemp

        # dssStx = pd.DataFrame(STXDSS, dss.Lines.AllNames(), ['A', 'B', 'C'])
        # dssSrx = pd.DataFrame(SRXDSS, dss.Lines.AllNames(), ['A', 'B', 'C'])
        # print("OpenDSS Loads after solving, from dss.CktElement.Powers():")
        # print(dss.CktElement.Powers())
        # return dssV, dssI, dssStx, dssSrx
        return VDSS


def solve_with_dss_constant_loads(dss_file):
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
        dss.Solution.Cleanup()

        print(
            f'\nOpendss Iterations: {dss.Solution.Iterations()}\tOpendss Convergence: {dss.Solution.Convergence()}')

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

            ph = np.asarray(dss.CktElement.BusNames()[
                            0].split('.')[1:], dtype='int')-1
            Imn = np.asarray(dss.CktElement.Currents())/Ibase
            Imn = Imn[0:int(len(Imn)/2)]
            Imn = Imn[0:5:2] + 1j*Imn[1:6:2]
            IDSS[k1, ph] = Imn
        dssI = pd.DataFrame(IDSS, dss.Lines.AllNames(), ['A', 'B', 'C'])

        STXDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')
        SRXDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')

        for k1 in range(len(dss.Lines.AllNames())):
            dss.Lines.Name(dss.Lines.AllNames()[k1])
            ph = np.asarray(dss.CktElement.BusNames()[
                            0].split('.')[1:], dtype='int')-1
            Sk = np.asarray(dss.CktElement.Powers())/(Sbase/1000)
            STXtemp = Sk[0:int(len(Sk)/2)]
            SRXtemp = Sk[int(len(Sk)/2):]
            STXtemp = STXtemp[0:5:2] + 1j*STXtemp[1:6:2]
            SRXtemp = -(SRXtemp[0:5:2] + 1j*SRXtemp[1:6:2])
            STXDSS[k1, ph] = STXtemp
            SRXDSS[k1, ph] = SRXtemp

        dssStx = pd.DataFrame(STXDSS, dss.Lines.AllNames(), ['A', 'B', 'C'])
        dssSrx = pd.DataFrame(SRXDSS, dss.Lines.AllNames(), ['A', 'B', 'C'])
        print("OpenDSS Loads after solving, from dss.CktElement.Powers():")
        print(dss.CktElement.Powers())
        return dssV, dssI, dssStx, dssSrx
