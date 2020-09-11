# # Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 10 September 2020
# Compare running time of opendss' powerflow solver to powerflowpy/fbs

import sys
import powerflowpy
from powerflowpy.fbs import *
import opendssdirect as dss
import time
import functools
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

def solver_timer(solver_func):
    """Decorator that returns the elapsed time when running solver_func n times"""
    @functools.wraps(solver_func)
    def wrapper(file, n):
        t1= time.perf_counter()
        for i in range(n):
            solver_func(file, n) # run once
        t2= time.perf_counter()
        return t2 - t1
    return wrapper

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

@solver_timer
def time_dss(dss_file, n):
    """wrap dss solver in a timer"""
    solve_with_dss(dss_file)


@solver_timer
def time_fbs(dss_file, n):
    """ wrap powerflowpy.fbs in a timer"""
    fbs(dss_file)

if __name__ == '__main__':
    """
    Compare opendss's solver with fbs and plot the results.
    Code based on Bernard Knasmueller's accurate-time-measurements-python at this github link: https://gitlab.com/bernhard.knasmueller/accurate-time-measurements-python/-/blob/master/main.py
    We use a straight-line fitting technique to eliminate systemic measurement errors
    The approach is described in this article: https://knasmueller.net/measure-code-execution-time-accurately-in-python
    And this paper: https://uwaterloo.ca/embedded-software-group/sites/ca.embedded-software-group/files/uploads/files/ieee-esl-precise-measurements.pdf
    """
    num_trials = 20
    dss_file = sys.argv[1]

    # initialize duration arrays, to avoid time lost to dynamic array resizing
    dss_durations = [-1] * num_trials
    fbs_durations = [-1] * num_trials

    for n in range(0,num_trials): #solve powerflow n...num_trials times with opendss
        dss_durations[n] = time_dss( dss_file, n)

    for n in range(0, num_trials):  # solve powerflow n...num_trials times with fbs
        fbs_durations[n] = time_fbs(dss_file, n)

    # change s to ms
    dss_durations = np.array([ 1000*s for s in dss_durations ])
    fbs_durations = np.array([ 1000*s for s in fbs_durations ])
    # initialize array for x axis, x = num trials
    x = np.array([i for i in range(1, num_trials + 1)])

    #fit data to a line
    dss_m, dss_b = np.polyfit(x, dss_durations, 1)
    fbs_m, fbs_b = np.polyfit(x, fbs_durations, 1)
    dss_f = lambda x: dss_m *x  + dss_b
    fbs_f = lambda x: fbs_m *x  + fbs_b

    # plot dss vs fbs
    fig, ax = plt.subplots()
    ax.set_xlabel("number of times solved")
    ax.set_ylabel("execution time(ms)")
    ax.plot(x, dss_durations, 'go', label = 'opendss')
    ax.plot(x, dss_f(x), color = 'green')
    ax.plot(x, fbs_durations, 'bo', label = 'fbs')
    ax.plot(x, fbs_f(x), color = 'blue')
    ax.legend()
    ax.grid(True)
    plt.savefig('dss_fbs_compare.png')
    plt.show()






