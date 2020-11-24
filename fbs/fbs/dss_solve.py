# # Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 10 September 2020
# run opendss' solver on a dss file. Copied form '20180601/opendss_nonvec_test_comparison.ipynb'

import opendssdirect as dss # type: ignore
import numpy as np # type: ignore
import pandas as pd # type: ignore
from typing import Tuple
from . utils import set_zip_values


def solve_with_dss(dss_file: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Set up dss object, call solve, and return the solution.
    """
    setup(dss_file)
    solve()
    return get_solution()


def setup(dss_file: str) -> None:
    """
    Initialize the network with a dss file, with one call to 'Redirect'.
    Set loads on the network.
    """
    dss.run_command('Redirect ' + dss_file)
    # originalSteps = dss.Solution.Number() # used to save original steps
    dss.Solution.Mode(1)
    dss.Solution.Number(1)
    dss.Solution.StepSize(1)
    dss.Solution.ControlMode(-1)
    dss.Solution.MaxControlIterations(1000000)
    dss.Solution.MaxIterations(30000)


def solve() -> None:
    """
    Call SolveSnap() to solve powerflow for one timestep.
    Save the loads before solving.
    """

    # Solve power flow with OpenDSS file, and updating loads before solving once
    # Code based on https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Version8/Distrib/Examples/Python/Python-to-OpenDSS%20Control%20Interface.pdf
    # save original loads

    # set_zip_values(dss)
    orig_loads_data = dss.utils.loads_to_dataframe()
    orig_loads_data = orig_loads_data.transpose()

    # this would run for a default value of originalSteps = 100
    # TODO: confirm with Shammya that this is correct, specifically that we do one timestep and call Solution.SolveSnap()
    for stepNumber in range(1):  # simulate timesteps

        # set loads
        for load_name in dss.Loads.AllNames():
            load_data = orig_loads_data[load_name]
            dss.Loads.Name(load_name)
            dss.Loads.kW(load_data['kW'])
            dss.Loads.kvar(load_data['kvar'])


        # run solve for this timestep
        dss.Solution.SolveSnap()
        dss.Solution.FinishTimeStep()

        if not dss.Solution.Converged():
            print('Initial Solution Not Converged. Check Model for Convergence')
        else:
            # Doing this solve command is required for GridPV, that is why the monitors
            # go under a reset process
            dss.Monitors.ResetAll()


def get_solution() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Input: A dss object with a solved network.
    Output: Network solution paramters V, I, Stx, Srx in pandas DataFrames.
    Also prints iterations and convergence to stdout∏.
    """
    print(f'\nOpendss Iterations: {dss.Solution.Iterations()}\tOpendss Convergence: {dss.Solution.Convergence()}')
    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0
    Ibase = Sbase/Vbase

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
        upstream_bus = dss.CktElement.BusNames()[0].split('.')[0]
        ph = np.asarray(dss.CktElement.BusNames()[0].split('.')[1:], dtype='int')-1
        # get upstream bus voltage
        dss.Circuit.SetActiveBus(upstream_bus)
        upstream_Vbase = dss.Bus.kVBase()*1000
        Ibase = Sbase/upstream_Vbase
        Imn = np.asarray(dss.CktElement.Currents())/Ibase
        Imn = Imn[0:int(len(Imn)/2)]
        Imn = Imn[0:5:2] + 1j*Imn[1:6:2]
        IDSS[k1, ph] = Imn
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

    loads = dict()
    for name in dss.Loads.AllNames():
        dss.Loads.Name(name)
        loads[name] = dss.CktElement.Powers()

    return dssV, dssI, dssStx, dssSrx, loads
