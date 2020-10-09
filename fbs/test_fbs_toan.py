from powerflowpy.utils import init_from_dss
from powerflowpy.fbs import fbs, get_solution as get_fbs_solution
from powerflowpy.dss_solve import solve_with_dss, get_solution as get_dss_solution

# Initialize the Network object from a dss file
# Loads will be initialized based on the dss file
dss_file_path = '/home/toanngo/Documents/GitHub/cigar-document/ceds-cigar/pycigar/data/ieee13busdata/ieee13.dss'
#dss_file_path = '/home/toanngo/Documents/GitHub/LinDist3Flow/nr3/nr3/data/06node_threephase_unbalanced.dss'
network = init_from_dss(dss_file_path)

# Set number of times to run fbs
timesteps = 10

# pick values for kW, kvar
kW = 1
kvar = 1 

for i in range(timesteps):
    fbs(network) # solve network once
    
    # update loads before next solve
    for load in network.get_loads():
        load.set_kW(kW)
        load.set_kvar(kvar)