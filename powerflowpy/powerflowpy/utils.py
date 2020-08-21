import numpy as np
import opendssdirect as dss
import sys
import scipy as sp
import scipy.sparse.csgraph as cg
import scipy.sparse.linalg as linalg
from typing import Iterable, List
from network import Node

def init_from_dss(network: Network, dss_fp: str) -> None:
    """define Network attributes from a dss file"""
    dss.run_command('Redirect ' + dss_fp)
    # set base values
    # make nodes
    for node_name in dss.Circuit.AllNodeNames():
        name, phase  = node_name.split('.')
        # get the node corresponding to this name, or make a new one
        if name not in network.nodes.keys():
            network.nodes[name] = Node(name)
        node = network.nodes[name]
        phase_idx = int(phase) - 1 # shift to 0-indexing
        self.phases[phase_idx] = 1 # add this phase to the node
    #make lines
    line_codes = dss.LineCodes.AllNames()
    lengths = [i() for i in dss.utils.Iterator(dss.Lines, 'Length')]
    line_zip = zip(line_codes, lengths)
    
    for line_code,length in line_zip:
        tx, rx = line_code.split('_')
        if (tx, rx) not in network.lines.keys():
            network.lines[(tx,rx)] = Line((tx,rx))
        line = network.lines[(tx,rx)]
        line.length = length
    # TODO: figure out how to get 3x3 impedance per unit matrix. 
    # probably something using dss.Lines, 'XMatrix' and 'RMatrix'
    

def is_acyclic(dss_fp: str) -> bool:
    #TODO Might be nice to replace with a cycle checker independent of opendssdirect
    dss.run_command('Redirect ' + dss_fp)
    return dss.Topology.NumLoops == 0
    


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
