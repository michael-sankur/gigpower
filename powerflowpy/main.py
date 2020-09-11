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

def solver_timer(solver_func, n):
    """Decorator that times running solver_func n times"""
    @functools.wraps(solver_func)
    def wrapper(*args, **kwargs):
        t1= time.perf_counter()
        for i in range(n):
            solver_func(*args, **kwargs) # run once
        t2= time.perf_counter()
        return t2 - t1
    return wrapper

@solver_timer
def time_dss(dss_file):
    """wrap dss solver in a timer"""
    dss.run_command('Redirect ' + dss_file)
    dss.Solution.Solve() # TODO: determine which command to use, this one or Solution.SolvePFlow()?

@solver_timer
def time_fbs(dss_file):
    """ wrap powerflowpy.fbs in a timer"""
    fbs(dss_file)

if __name__ == '__main__':
    """
    Compare opendss's solver with fbs and plot the results.
    We use a straight-line fitting technique to eliminate systemic measurement errors
    The approach is described in this article: https://knasmueller.net/measure-code-execution-time-accurately-in-python
    And this paper: https://uwaterloo.ca/embedded-software-group/sites/ca.embedded-software-group/files/uploads/files/ieee-esl-precise-measurements.pdf
    """
    num_trials = 200
    dss_file = sys.argv[1]

    # get the first datapoint for each funciton
    dss_durations = [ time_dss( dss_file, 1) ]
    dss_durations = [ time_fbs( dss_file, 1) ]

    for n in range(2,num_trials+1): #solve powerflow n...num_trials times with opendss
        dss_durations.append( time_dss(dss_file,  n) )

    for trial in num_trials(200):  # solve powerflow num_trials times with fbs
        total_fbs_time += time_dss(dss_file)

    # return the estimated total time



