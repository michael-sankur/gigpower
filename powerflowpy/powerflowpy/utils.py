import numpy as np
import opendssdirect as dss
import sys
from typing import Iterable, List, Any
from . network import Network, Node, Line, Load, Controller, Capacitor
import scipy.linalg as spla

def init_from_dss(dss_fp: str) -> None:
    """define a Network attributes from a dss file"""
    dss.run_command('Redirect ' + dss_fp)
    network = Network()

    # set base values
    dss.Solution.Solve()
    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0
    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase

    # make Nodes
    for node_name in dss.Circuit.AllNodeNames():
        name, phase  = node_name.split('.')
        # get the node corresponding to this name, or make a new one
        if name not in network.nodes.keys():
            network.nodes[name] = Node(name)
        node = network.nodes[name]
        phase_idx = int(phase) - 1 # shift to 0-indexing
        node.phases[phase_idx] = 1 # add this phase to the node
        # initalize this node's adjacency list to the empty list
        # note: this means that every node has an entry. Nodes with no children will hav an empty list.
        network.adj[name] = []

    #make Lines
    all_lines_data = dss.utils.lines_to_dataframe().transpose() # get dss line data indexed by line_code
    line_codes = all_lines_data['LineCode'].keys()

    for line_code in line_codes:
        line_data = all_lines_data[line_code]
        tx, *tx_phases = line_data['Bus1'].split('.')
        rx, *rx_phases = line_data['Bus2'].split('.')
        if tx_phases != rx_phases:
            raise ValueError(f'Tx phases do not match Rx phases for line {line_code}')
        line = Line((tx, rx)) # initialize line
        for phase in tx_phases: # set phases according to tx
            line.phases[get_phase_idx(phase)] = 1

        network.lines[(tx,rx)] = line # add line to network.line
        # add directed line to adjacency list, adj[tx] += rx
        network.adj[tx].append(rx)

        #parse line attributes from dss line data
        line.name = line_code
        line.length = line_data['Length']
        line.FZpu = get_Z_from_Y(line_data['YPrim'], line.phases)


    # make Loads
    load_names = dss.Loads.AllNames()
    for load_name in load_names:
        # TODO: handle multiple loads on a node
        node_name, phase_char, load_idx = load_name.split('_')[1:]
        try:
            node = network.nodes[node_name]
        except KeyError:
            print("Node assigned to load not defined for this network.")
        load = Load(load_name + '_' + load_idx)
        load.phases[get_phase_idx(phase_char)] = 1 # indicate this load's phase
        node.load = load #assign this load to its node
    #TODO: get aPQ, aI, ppu, qpu, spu for each load
    #TODO: figure out how to get connection information (wye or delta).
    # Is this info determined by Capacitors?
    # Or can we use Loads.IsDelta?

    # make Controllers
    # TODO: implement this. No idea how opendssdirect maps this info. Which class is it even?

    # make Capacitors
    cap_names = dss.Capacitors.AllNames()
    for cap_name in cap_names:
        cap = Capacitor(cap_name)
        #TODO: figure out how to parse these names from an example
    return network

def get_phase_idx(phase_char: str) -> int:
    """
    helper function to turn a phase letter into an index, where 'a' = 0
    """
    if phase_char in ['a', 'b', 'c']:
        return ord(phase_char.lower()) - ord('a')
    elif phase_char in ['1', '2', '3']:
        return int(phase_char) - 1
    else:
        raise ValueError(f'Invalid argument for get_phase_idx {phase_char}')

def get_Z_from_Y(yp: Iterable, phases : List ) -> Iterable:
    """
    helper function to get the Z matrix from the flattened Yprim matrix.
    Returns a numpy array.
    """
    #TODO: use phases to pad Y matrix as necessary [phases = [1 1 0]]
    yp = np.asarray(yp, dtype='complex')
    yp = yp[0:-1:2] + 1j*yp[1::2]
    yp = np.reshape(yp, (int(yp.shape[0]**(1/2)), int(yp.shape[0]**(1/2))))
    return spla.inv(yp)
