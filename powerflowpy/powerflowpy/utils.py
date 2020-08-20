import numpy as np
import opendssdirect as dss
import sys
import scipy as sp
import scipy.sparse.csgraph as cg
import scipy.sparse.linalg as linalg
from typing import Iterable, List


def is_acyclic(dss_fp: str) -> bool:
    #TODO Might be nice to replace with a cycle checker independent of opendssdirect
    dss.run_command('Redirect ' + dss_fp)
    return dss.Topology.NumLoops == 0
    
def get_network(dss_fp: str) -> Iterable :
    """ returns a Compressed Sparse Column (CSC) representation of the network"""
    dss.run_command('Redirect ' + dss_fp)
    return sp.sparse.csc_matrix(dss.YMatrix.getYsparse())

def get_node_names(dss_fp: str) -> List :
    """ 
    returns a 1x array of node_names corresponding to column indices of the CSC
    node_names[column_index] = 'name of node represented by column_index'
    """
    dss.run_command('Redirect ' + dss_fp)
    return dss.Circuit.AllNodeNames())

def get_bfs_order(network: Iterable) -> List:
    """ returns a 1x array of nodes in breadth-first-search order from the node at index 0 """
    #TODO: determind if the source will always be the first bus, or if we need to 
    #derive the source from the graph topology 

def find_source(network: Iterable) -> int:
    #TODO: It might be useful to define the root using the network topology.
    # Currently relying on opendss to identify the root as the sourcebus 
    pass