import sys
import powerflowpy
from powerflowpy.utils import init_from_dss
from powerflowpy.fbs import *
import opendssdirect as dss

if __name__ == '__main__':
    # Runs from command line with path to dss file as an argument
    dss_file = sys.argv[1]
    dss.run_command('Redirect ' + dss_file)
    dss.Solution.SolvePFlow()
    print(dss.Solution.TotalIterations())
    network = init_from_dss(dss_file)
    solution = fbs(dss_file)
    solution.print_solution()
