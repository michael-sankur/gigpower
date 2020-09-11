# # Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 10 September 2020
# run opendss' solver on a dss file

import opendssdirect as dss

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


