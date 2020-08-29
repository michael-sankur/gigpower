import numpy as np
import opendssdirect as dss
import sys
from typing import Iterable, List, Tuple, Any
from . network import Network, Node, Line, Load, Controller, Capacitor

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
            network.nodes[name].phases = [] # replace default tuple with empty list so we can mutate it
        node = network.nodes[name]
        # store phases as characters in phase list for now
        node.phases.append(phase)
        # Add node to adjacency lisst. note: this means that every node has an entry. Nodes with no children will have an empty list.
        network.adj[name] = []
    # iterate through nodes to parse phase lists
    for node in network.get_nodes():
        node.phases = parse_phases(node.phases)

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
        # set rx's parent to tx
        network.nodes.get(rx).parent = network.nodes.get(tx)

        #parse line attributes from dss line data
        line.name = line_code
        line.length = line_data['Length']
        line.FZpu = get_Z(line_data, line.phases)


    # make Loads
    all_loads_data = dss.utils.loads_to_dataframe().transpose()
    load_names = all_loads_data.keys()
    for load_name in load_names:
        # TODO: handle multiple loads on a node
        load_data = all_loads_data[load_name]
        node_name, phase_char, load_idx = load_name.split('_')[1:]
        try:
            node = network.nodes[node_name]
        except KeyError:
            print(f"This load's node has not been defined. Load name: {load_name}, Node name: {node_name}")
        load = Load(load_name + '_' + load_idx)
        load.phases = parse_phases([phase_char])
        node.load = load #assign this load to its node
        load.ppu = load_data['kW'] / 1000
        load.qpu = load_data['kvar']  / 1000
        load.spu = load.ppu + 1j*load.qpu
        load.conn =   'delta' if load_data['IsDelta'] else 'wye' #TODO: figure out if delta/wye are mutually exclusive
        #TODO: set load.type
        #TODO: get aPQ, aI for each load

    # make Controllers
    # TODO: implement this. No idea how opendssdirect maps this info. Which class is it even?

    # make Capacitors
    all_cap_data = dss.utils.capacitors_to_dataframe().transpose()
    cap_names = all_cap_data.keys()
    for cap_name in cap_names:
        cap_data = all_cap_data[cap_name]
        cap = Capacitor(cap_name)
        # TODO: figure out if delta/wye are mutually exclusive]
        cap.conn = 'delta' if cap_data['IsDelta'] else 'wye'
        cap.cappu = cap_data['kvar'] * 1000 / self.Sbase
        #TODO: get phases on capacitor

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

def get_Z(dss_data: Any, phase_list : Tuple ) -> Iterable:
    """
    helper function to get the Z matrix from dss.lines.to_dataframe()
    Returns an ndarray.
    """
    #TODO: how to make this per unit length
    num_phases = phase_list.count(True)
    RM = np.asarray(dss_data['RMatrix'])
    XM = np.asarray(dss_data['XMatrix'])
    # reshape based on phases
    phases = int(dss_data['Phases'])
    ZM = RM + 1j*XM
    ZM = np.reshape(ZM, (ZM.shape[0]//phases, phases))  # reshape
    # pad the Z matrix
    z_padded = np.zeros((3, 3), dtype=complex)
    z_vals = iter(ZM.flatten())
    for row_idx in range(3):
        for col_idx in range(3):
            if phase_list[row_idx] and phase_list[col_idx]:
                try:
                    z_padded[row_idx][col_idx] = next(z_vals)
                except StopIteration:
                    (f"There is a mismatch in phases between line {line_name} and the Z matrix")
    return z_padded

def mask_phases(matrix: Iterable, phases: tuple) -> Iterable:
    """
    Zeroes out values in input matrix for phases set to FALSE in the phases tuple.
    Input:
        matrix: a 3x3 ndarray
        phases: a tuple of booleans corresponding to phases to set on this matrix (A: T/F, B: T/F, C: T/F)
    Output:
        input matrix with 0's for all row/column indices corresponding to phases set to FALSE
    """
    # create a 3x3 phase matrix of 1's and 0's base on phases
    phase_matrix = np.zeros((3, 3), dtype=complex)
    for row_idx in range(3):
        for col_idx in range(3):
            if phases[row_idx] and phases[col_idx]:
                phase_matrix[row_idx, col_idx] = 1

    return np.matmul(matrix, phase_matrix)
