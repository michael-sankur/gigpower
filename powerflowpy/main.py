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
    print("Opendss Solution")
    print(f"iterations: {dss.Solution.TotalIterations()}")
    nodes = dss.Circuit.AllBusVolts()
    lines = dss.Circuit.YCurrents()
    print(nodes)
    print(lines)

    network = init_from_dss(dss_file)
    solution = fbs(dss_file)
    print("\nPython FBS Solution")
    solution.print_solution()
    print()