# # Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 10 September 2020
# Compare running time of opendss' powerflow solver to fbs/fbs

import click
from pathlib import Path
from fbs.fbs import fbs
from fbs.utils import init_from_dss
from fbs.dss_solve import setup, solve
from fbs.network import Network
import time
import functools
import numpy as np # type: ignore
import csv
import matplotlib.pyplot as plt # type: ignore
from typing import Any
plt.style.use('seaborn-whitegrid')


def solver_timer(solver_func : Any) -> Any:
    """Decorator that returns the elapsed time when running solver_func n times"""
    @functools.wraps(solver_func)
    def wrapper(network : Network, n : int) -> Any:
        t1 = time.perf_counter()
        for i in range(n):
            solver_func(network, i)  # run once
        t2 = time.perf_counter()
        return t2 - t1
    return wrapper

@solver_timer
def time_dss(dss_file: str, i: int) -> None:
    """wrap dss solver in a timer"""
    setup(dss_file)
    solve()
    return None


@solver_timer
def time_fbs(network: Network, i: int) -> None :
    """ wrap powerflowpy.fbs in a timer"""
    fbs(network)
    return None

# set up command line arguments to this script
@click.command( help='Compare running times between opendss and fbs. Please provide a full or relative path to a dss file.')
@click.argument('dss_file')
@click.option('--output', metavar='PATH', help='write output files to this directory')
@click.option('--trials', default = '100', help = 'Specify max number of trials to run. Default is 20.')

def main(dss_file: str, trials: int, output: str) -> None:
    """
    Compare opendss's solver with fbs and plot the results.
    Code based on Bernard Knasmueller's accurate-time-measurements-python at this github link: https://gitlab.com/bernhard.knasmueller/accurate-time-measurements-python/-/blob/master/main.py
    We use a straight-line fitting technique to eliminate systemic measurement errors
    The approach is described in this article: https://knasmueller.net/measure-code-execution-time-accurately-in-python
    And this paper: https://uwaterloo.ca/embedded-software-group/sites/ca.embedded-software-group/files/uploads/files/ieee-esl-precise-measurements.pdf

    Sample command with default number of trials, and an output directory specified
    $ python3 fbs/main.py 'fbs/tests/06n3ph_rad_unbal/06node_threephase_radial_unbalanced.dss' --output ./fbs/tests/06n3ph_rad_unbal
    """
    num_trials = int(trials)

    # initialize duration arrays, to avoid time lost to dynamic array resizing
    dss_durations = [-1] * num_trials
    fbs_durations = [-1] * num_trials

    for n in range(0,num_trials): #solve powerflow n...num_trials times with opendss
        dss_durations[n] = time_dss( dss_file, n+1)

    # create a network objet for fbs to read
    network = init_from_dss(dss_file)
    for n in range(0, num_trials):  # solve powerflow n...num_trials times with fbs
        fbs_durations[n] = time_fbs(network, n+1)
    # change s to ms
    dss_durations = np.array([ 1000*s for s in dss_durations ])
    fbs_durations = np.array([ 1000*s for s in fbs_durations ])
    # initialize array for x axis, x = num trials
    x = np.array([i for i in range(1, num_trials + 1)])

    # if output is specified, write csv file
    if output:
        parent = Path(output)
        parent.mkdir(exist_ok = True)
        with parent.joinpath('fbs_v_dss.csv').open('w') as out:
            writer = csv.writer(out)
            writer.writerow([f'Input file: {dss_file}'])  # write path to input file
            writer.writerow(['n', 'fbs', 'dss'])  # write header row
            for row in zip( x, fbs_durations, dss_durations): # write timing data
                writer.writerow(row)
    # otherwise, write to stdout
    else:
        click.echo("n, fbs, dss")
        click.echo(zip(x, fbs_durations, dss_durations))

    #fit data to a line
    dss_m, dss_b = np.polyfit(x, dss_durations, 1)
    fbs_m, fbs_b = np.polyfit(x, fbs_durations, 1)
    dss_f = lambda x: dss_m *x  + dss_b
    fbs_f = lambda x: fbs_m *x  + fbs_b

    # plot dss vs fbs
    fig, ax = plt.subplots()
    plt.title(f'{dss_file}')
    ax.set_xlabel("number of times solved")
    ax.set_ylabel("execution time(ms)")
    ax.plot(x, dss_durations, 'go', label = 'opendss')
    ax.plot(x, dss_f(x), color='green', # type: ignore
            label=f"dss time/trial ~= {dss_m :.2f}ms")
    ax.plot(x, fbs_durations, 'bo', label = 'fbs')
    ax.plot(x, fbs_f(x), color='blue', label=f'fbs time/trial ~={fbs_m: .2f}ms') # type: ignore

    plt.xticks(np.arange(1, num_trials + 1)) # set x ticks to whole numbers
    ax.legend()
    ax.grid(True)
    if output:
        parent = Path(output)
        parent.mkdir(exist_ok=True)
        plt.savefig(parent.joinpath('dss_fbs_compare.png'))
    plt.show()
    print(f'FBS m: {fbs_m}ms')
    print(f'dss m: {dss_m}ms')

if __name__ == '__main__':
    main()


