# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Implement FBS to solve Power Flow for a radial distribution network. 

import sys
from typing import Iterable
from utils import is_acyclic, get_network, get_node_names

def fbs(dss_fp) -> None:

    if not is_acyclic(dss_fp):
        raise ValueError('Not a radial network.')

    network = get_network(dss_fp) # get network as a compressed matrix
    node_names = get_node_names(dss_fp) # get node names corresponding to network indices
    
    # Settheinitialvoltagetobebalancedthree-phasevoltage 
    # atthebuswhereislocatedatlineterminalorabranch intersection. 
    # Start the proposed procedure from the farthestbusobtainedbystep2.

if __name__ == '__main__':
    # run fbs from command line with path to dss file as an argument
    dss_file = sys.argv[1]
    # run fbs
    fbs( get_network(dss_file) )
