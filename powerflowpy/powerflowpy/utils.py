import numpy as np
import opendssdirect as dss
import sys
from typing import Iterable, List, Tuple, Any
from . network import Network, Node, Line, Load, Controller, Capacitor
import scipy.linalg as spla

def init_from_dss(dss_fp: str) -> None:
    """define a Network attributes from a dss file"""
    dss.run_command('Redirect ' + dss_fp)
    network = Network()

    # set base values
    # TODO: make it possible to set the base from a given bus
    dss.Solution.Solve()
    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0
    Ibase = Sbase/Vbase
    network.Zbase = Vbase/Ibase
    network.Vbase, network.Sbase, network.Ibase = Vbase, Sbase, Ibase

    # make Nodes
    for node_name in dss.Circuit.AllNodeNames():
        name, phase  = node_name.split('.')
        # get the node corresponding to this name, or make a new one
        if name not in network.nodes.keys():
            network.nodes[name] = Node(name)
        node = network.nodes[name]
        node.phases = parse_phases([phase])
        # note: this means that every node has an entry. Nodes with no children will hav an empty list.
        network.adj[name] = []

    #make Lines
    all_lines_data = dss.utils.lines_to_dataframe().transpose() # get dss line data indexed by line_code
    line_codes = all_lines_data.keys()

    for line_code in line_codes:
        line_data = all_lines_data[line_code]
        tx, *tx_phases = line_data['Bus1'].split('.')
        rx, *rx_phases = line_data['Bus2'].split('.')
        if tx_phases != rx_phases:
            raise ValueError(f'Tx phases do not match Rx phases for line {line_code}')
        line = Line((tx, rx), line_code) # initialize line
        for phase in tx_phases: # set phases according to tx
            line.phases = parse_phases(tx_phases)
        network.lines[(tx,rx)] = line # add line to network.line
        # add directed line to adjacency list, adj[tx] += rx
        network.adj[tx].append(rx)

        #parse line attributes from dss line data
        line.name = line_code
        line.length = line_data['Length']
        line.FZpu = get_Z_from_Y(line_data['Yprim'], line.phases)


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
        load.phases = parse_phases([phase_char])
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


def parse_phases(phase_char_lst):
    """
    helper function to return a list of phase characters into a boolean triple
    ex: ['1', '3'] -> (True, False True)
    ex: ['a','b'] -> (True, True, False)
    """
    phase_list = [False, False, False]
    for p in phase_char_lst:
        phase_list[get_phase_idx(p)] = True
    return tuple(phase_list)

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

def get_Z_from_Y(YP: Iterable, phase_list : Tuple ) -> Iterable:
    """
    helper function to get the Z matrix by inverting corresponding elements of the
    Yprim matrix available from dss.lines.to_dataframe()
    Returns an ndarray.
    """
    num_phases = phase_list.count(True)
    YP = np.asarray(YP)
    YP = YP[0:-1:2] + 1j*YP[1::2]  # reshape into pairs, real and complex
    # reshape based on phases
    YP = np.reshape(YP, (YP.shape[0]//num_phases, num_phases))
    YP = YP[0:num_phases*2:2]  # take the rows corresponding to ZPrime phases
    invY = spla.inv(YP)
    # pad the Y matrix
    yp_padded = np.zeros((3, 3), dtype=complex)
    yp_vals = iter(invY.flatten())
    for row_idx in range(3):
        for col_idx in range(3):
            if phase_list[row_idx] and phase_list[col_idx]:
                try:
                    yp_padded[row_idx][col_idx] = next(yp_vals)
                except StopIteration:
                    (f"There is a mismatch in phases between line {line_name} and the dss.YPrim matrix. \n Here is the line data: {this_line}")
    return yp_padded
