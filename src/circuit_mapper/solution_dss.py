# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: March 23 2021
# Create SolutionDSS class, a Solution that uses opendss to solve. 

import numpy as np  # type: ignore
from typing import List, Dict, Iterable
from . solution import Solution
from . circuit import Circuit
from . utils import set_zip_values_dss, parse_phases, pad_phases
import opendssdirect as dss  # type: ignore
import numpy as np
import pandas as pd
from typing import Tuple
from collections import defaultdict


class SolutionDSS(Solution):

    def __init__(self, dss_fp: str):
        super().__init__(dss_fp) # saves a Circuit object
        self.dss_fp = dss_fp

    def solve(self):
        """
        Set up dss object, call solve_interface, and return the solution.
        """
        self.setup()
        self.solve_interface()
        self.save_solution()

    def setup(self):
        """
        Initialize the network with a dss file, with one call to 'Redirect'.
        Set loads on the network.
        """
        dss.run_command('Redirect ' + self.dss_fp)
        # originalSteps = dss.Solution.Number() # used to save original steps
        dss.Solution.Mode(1)
        dss.Solution.Number(1)
        dss.Solution.StepSize(1)
        dss.Solution.ControlMode(-1)
        dss.Solution.MaxControlIterations(1000000)
        dss.Solution.MaxIterations(30000)
        self.dss = dss

    def solve_interface(self):
        """
        Call SolveSnap() to solve powerflow for one timestep.
        Save the loads before solving.
        """
        # Solve power flow with OpenDSS file, and updating loads before solving once
        # Code based on https://sourceforge.net/p/electricdss/code/HEAD/tree/
        # trunk/Version8/Distrib/Examples/Python/Python-to-OpenDSS%20Control%20Interface.pdf
        
        # save original loads
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
            SlackBusVoltage = 1.00
            dss.Vsources.PU(SlackBusVoltage)

            set_zip_values_dss(dss, Circuit.ZIP_V)
            # run solve for this timestep
            dss.Solution.SolveSnap()
            dss.Solution.FinishTimeStep()

            if not dss.Solution.Converged():
                print('Initial Solution Not Converged. Check Model for Convergence')
            else:
                # Doing this solve command is required for GridPV, that is why the monitors
                # go under a reset process
                dss.Monitors.ResetAll()
        
        # save convergance and iteracitons
        self.iterations = self.dss.Solution.Iterations()
        self.convergence_diff = self.dss.Solution.Convergence()

    def save_solution(self):
        """
        Saves V, I, STX, SRX, and sV from current state of dss
        """
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
        self.V = VDSS
        self.Vmag = np.abs(VDSS)

        IDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')
        for k1 in range(len(dss.Lines.AllNames())):
            dss.Lines.Name(dss.Lines.AllNames()[k1])
            upstream_bus = dss.CktElement.BusNames()[0].split('.')[0]
            ph = np.asarray(dss.CktElement.BusNames()[0].split('.')[1:], dtype='int')-1
            if len(ph) == 0:  # if the name doesn't have phases, assume all three phases
                ph = [0, 1, 2]
            # get upstream bus voltage
            dss.Circuit.SetActiveBus(upstream_bus)
            upstream_Vbase = dss.Bus.kVBase()*1000
            Ibase = Sbase/upstream_Vbase
            Imn = np.asarray(dss.CktElement.Currents())/Ibase
            Imn = Imn[0:int(len(Imn)/2)]
            Imn = Imn[0:5:2] + 1j*Imn[1:6:2]
            IDSS[k1, ph] = Imn
        self.I = IDSS

        STXDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')
        SRXDSS = np.zeros((len(dss.Lines.AllNames()), 3), dtype='complex')

        for k1 in range(len(dss.Lines.AllNames())):
            dss.Lines.Name(dss.Lines.AllNames()[k1])
            ph = np.asarray(dss.CktElement.BusNames()[0].split('.')[1:], dtype='int')-1
            if len(ph) == 0:  # if the name doesn't have phases, assume all three phases
                ph = [0, 1, 2]
            Sk = np.asarray(dss.CktElement.Powers())/(Sbase/1000)
            STXtemp = Sk[0:int(len(Sk)/2)]
            SRXtemp = Sk[int(len(Sk)/2):]
            STXtemp = STXtemp[0:5:2] + 1j*STXtemp[1:6:2]
            SRXtemp = -(SRXtemp[0:5:2] + 1j*SRXtemp[1:6:2])
            STXDSS[k1, ph] = STXtemp
            SRXDSS[k1, ph] = SRXtemp

        self.Stx = STXDSS
        self.Srx = SRXDSS
        self.sV = self.get_bus_powers().values / 1000

    def get_V_mag(self) -> pd.DataFrame:
        """Return VMag per bus in a DataFrame"""
        buses = self.dss.Circuit.AllBusNames()
        data = np.zeros((len(buses), 3), dtype=float)

        for i, b in enumerate(buses):
            self.dss.Circuit.SetActiveBus(b)
            ph = np.asarray(self.dss.Bus.Nodes(), dtype='int')-1
            Vtemp = np.asarray(self.dss.Bus.puVmagAngle())[0::2]
            data[i, ph] = Vtemp

        return pd.DataFrame(data, buses, ['A', 'B', 'C'])


    def get_nominal_bus_powers(self) -> pd.DataFrame:
        """
        Return total nominal power for loads and capacitors, summed by bus
        Capacitors' signs must be reversed from dss.Capacitors.kvar(); we must
        SUBTRACT capacitors
        """
        dss = self.dss
        bus_names = dss.Circuit.AllBusNames()
        bus_idx_dict = { b: i for i, b in enumerate(bus_names)}
        data = np.zeros((len(bus_names), 3), dtype=complex)
        
        # add all load powers to load data, based on bus index
        for load in dss.Loads.AllNames():
            dss.Loads.Name(load)
            bus_name = dss.CktElement.BusNames()[0]
            bus_name, bus_phase = bus_name.split('.')[0], bus_name.split('.')[1:]
            if len(bus_phase) == 0:
                bus_phase.extend(['1', '2', '3'])
            bus_idx = bus_idx_dict[bus_name]
            phases = parse_phases(list(bus_phase))
            dist_real = dss.Loads.kW() / len(bus_phase)
            dist_imag = dss.Loads.kvar() / len(bus_phase)
            data[bus_idx, phases] += (dist_real + 1j * dist_imag)

        # add all capacitor powers to load data, based on bus index
        for cap in dss.Capacitors.AllNames():
            dss.Capacitors.Name(cap)
            bus_name = dss.CktElement.BusNames()[0]
            bus_name, bus_phase = bus_name.split('.')[0], bus_name.split('.')[1:]
            if len(bus_phase) == 0:
                bus_phase.extend(['1', '2', '3'])
            bus_idx = bus_idx_dict[bus_name]
            phases = parse_phases(list(bus_phase))
            real, imag = 0, dss.Capacitors.kvar() / phases.count(True)
            data[bus_idx, phases] += (real - 1j * imag)

        return pd.DataFrame(data, bus_names, ['A', 'B', 'C'])


    def get_load_powers(self) -> pd.DataFrame:
        """
        """
        dss = self.dss
        load_names = dss.Loads.AllNames()
        cols = ['A', 'B', 'C']
        load_powers_by_bus = {name: np.zeros(3, dtype=complex) for name in dss.Circuit.AllBusNames()}

        for name in load_names:
            dss.Loads.Name(name)
            bus_name = dss.CktElement.BusNames()[0]
            bus_name, bus_phase = bus_name.split('.')[0], bus_name.split('.')[1:]
            if len(bus_phase) == 0:
                bus_phase.extend(['1', '2', '3'])
            phases = parse_phases(list(bus_phase))
            sparse = np.asarray(dss.CktElement.Powers())
            # pull out real and imaginary components, pad by phase, ignore neutral
            real = pad_phases(sparse[0:5:2], (3, ), phases)
            imag = pad_phases(sparse[1:6:2], (3, ), phases)
            load_powers_by_bus[bus_name] += real + 1j * imag

        return pd.DataFrame.from_dict(load_powers_by_bus, orient='index', dtype=complex, columns=cols)


    def get_capacitor_powers(self) -> pd.DataFrame:
        """
        Query opendss for Capacitor powers by bus
        helper method to getTotalNodePowers_CktElement
        """
        dss = self.dss
        cap_names = dss.Capacitors.AllNames()
        cols = ['A', 'B', 'C']
        cap_powers_by_bus = {b: np.zeros(3, dtype=complex) for b in dss.Circuit.AllBusNames()}

        # add all capacitors to load data, based on bus index
        for name in cap_names:
            dss.Capacitors.Name(name)
            bus_name = dss.CktElement.BusNames()[0]
            bus_name, bus_phase = bus_name.split('.')[0], bus_name.split('.')[1:]
            if len(bus_phase) == 0:
                bus_phase.extend(['1', '2', '3'])
            phases = parse_phases(list(bus_phase))
            sparse = np.asarray(dss.CktElement.Powers())
            # pull out real and imaginary components, pad by phase, ignore neutral
            real = pad_phases(sparse[0:5:2], (3, ), phases)
            imag = pad_phases(sparse[1:6:2], (3, ), phases)
            cap_powers_by_bus[bus_name] += real + 1j * imag

        return pd.DataFrame.from_dict(cap_powers_by_bus, orient='index', dtype=complex, columns=cols)

    def get_data_frame(self, param: str) -> pd.DataFrame:
        """
        overrides super() to handle Lines 
        """
        try:
            element_group, cols, data_type = self.__class__.SOLUTION_PARAMS.get(param)
            if element_group == 'lines':
                index = self.dss.Lines.AllNames()
            else:
                index = getattr(self.circuit, element_group).all_names()
            data = getattr(self, param)
            return pd.DataFrame(data=data, index=index, columns=cols, dtype=data_type)
        except KeyError:
            print(f"Not a valid solution parameter. Valid parameters: \
                  {self.__class__.SOLUTION_PARAMS.keys()}")