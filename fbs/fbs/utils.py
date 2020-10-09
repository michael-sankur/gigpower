import numpy as np # type: ignore
import opendssdirect as dss # type: ignore
from typing import Iterable, List, Any
from . network import Network, Node, Line, Load, Capacitor

def init_from_dss(dss_fp: str) -> Network:
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
        node.phases = parse_phases(node.phases) #type: ignore

    #make Lines
    all_lines_data = dss.utils.lines_to_dataframe().transpose() # get dss line data indexed by line_code
    line_codes = all_lines_data.keys()

    for line_code in line_codes:
        line_data = all_lines_data[line_code]

        # Handle the typical case where bus names tell you the phases, e.g. 'Bus1.1.2.'
        if "." in line_data['Bus1']:
            tx, *tx_phases = line_data['Bus1'].split('.')
            rx, *rx_phases = line_data['Bus2'].split('.')          
            if tx_phases != rx_phases:
                raise ValueError(f'Tx phases do not match Rx phases for line {line_code}')
            
        else:  # Otherwise, if all 3 phases are present, define the line for all 3 phases.
            tx, rx = line_data['Bus1'], line_data['Bus2']
            tx_phases = ['1', '2', '3']

        line = Line((tx, rx), line_code)  # initialize line
        # set phases according to tx
        line.phases = parse_phases(tx_phases)
        network.lines[(tx, rx)] = line  # add line to network.line
        # add directed line to adjacency list, adj[tx] += rx
        network.adj[tx].append(rx)
        tx_node, rx_node = network.nodes.get(tx), network.nodes.get(rx)
        # Make sure that rx does not yet have a parent assigned. If so, set rx's parent to tx
        if rx_node.parent:  # type: ignore
            raise ValueError(f"Not a radial network. Node {rx_node.name} has more than one parent.")  # type: ignore
        else:
            rx_node.parent = tx_node  # type: ignore

        # parse line attributes from dss line data
        line.name = line_code
        line.length = line_data['Length']
        fz_mult = 1 / network.Zbase * line.length

        line.FZpu = get_Z(line_data, line.phases, fz_mult)

    # make Loads
    all_loads_data = dss.utils.loads_to_dataframe().transpose()
    load_names = all_loads_data.keys()
    for load_name in load_names:
        load_data = all_loads_data[load_name]
        #node_name, phase_chars, load_idx = load_name.split('_')[1:]
        dss.Loads.Name(load_name)
        bus_phase = dss.CktElement.BusNames()[0].split('.')
        if len(bus_phase) == 1:
            bus_phase.extend(['1', '2', '3'])
        phase_chars = ''
        for i in range(dss.CktElement.NumPhases()):
            if bus_phase[i+1] == '1':
                phase_chars += 'a'
            if bus_phase[i+1] == '2':
                phase_chars += 'b'
            if bus_phase[i+1] == '3':
                phase_chars += 'c'
        node_name = bus_phase[0]
        try:
            node = network.nodes[node_name]
        except KeyError:
            print(f"This load's node has not been defined. Load name: {load_name}, Node name: {node_name}")
        load = Load(load_name)
        load.phases = parse_phases(list(phase_chars))
        # save kw and kvar
        load.kW = load_data['kW']
        load.kvar = load_data['kvar']

        # divide ppu and qpu by number of phases
        ppu = load.kW / 1000/ load.phases.count(True)
        qpu = load.kvar / 1000 / load.phases.count(True)

        # set aPQ, aI, aZ
        # TODO: get aPQ, aI for each load
        load.aPQ = np.ones(3)
        load.aI = np.zeros(3)
        load.aZ = np.zeros(3)
        # set load's ppu, qpu, and spu as a 3x1 based on load's phases
        load.ppu = np.asarray( [ppu if phase else 0 for phase in load.phases])
        load.qpu = np.asarray( [qpu if phase else 0 for phase in load.phases])
        load.spu = load.ppu + 1j*load.qpu
        load.conn =   'delta' if load_data['IsDelta'] else 'wye'
        # add a pointer to this load to the network
        network.loads[load_name] = load
        load.node = node # assign the node to the load
        node.loads.append(load)  # add this load to its node's load list
        #TODO: set load.type
    # sum all loads on each node and store on each node, to avoid re-calculating during
    # fbs.update-voltage-dependent-load()
    for node in network.nodes.values():
        for load in node.loads:
            node.sum_spu = np.add(node.sum_spu, load.spu)

    # make Controllers
    # TODO: implement this. No idea how opendssdirect maps this info. Which class is it even?

    # make Capacitors
    all_cap_data = dss.utils.capacitors_to_dataframe().transpose()
    cap_names = all_cap_data.keys()
    for cap_name in cap_names:
        if cap_name != '':
            cap_data = all_cap_data[cap_name]
            node_name, phase_chars, cap_idx = cap_name.split('_')[1:]
            try:
                node = network.nodes[node_name]
            except KeyError:
                print(
                    f"This cap's node has not been defined. Cap name: {cap_name}, Node name: {node_name}")
            cap = Capacitor(cap_name)
            cap.phases = parse_phases(list(phase_chars))
            cap.conn = 'delta' if cap_data['IsDelta'] else 'wye'
            cappu = cap_data['kvar'] * 1000 / network.Sbase / len(cap.phases) # TODO: confirm that cappu is divided by num phases
            cap.cappu = np.asarray([cappu if phase else 0 for phase in cap.phases])
            network.capacitors[ cap_name ] = cap # add a pointer to this cap to the network
            node.capacitors.append(cap) # add capacitor to it's node's cap list
    # sum all cappu on each node and store on each node, to avoid re-calculating during
    # fbs.update-voltage-dependent-load()
    for node in network.nodes.values():
        for cap in node.capacitors:
            node.sum_cappu = np.add(node.sum_cappu, cap.cappu)

    # # make Transformers
    # all_transformer_data = dss.utils.transformers_to_dataframe().transpose()
    # transformer_names = all_transformer_data.keys()
    # print(all_transformer_data)

    return network


def parse_phases(phase_char_lst : List[str]) -> List[bool]:
    """
    helper function to return a list of phase characters into a boolean triple
    ex: ['1', '3'] -> [True, False True]
    ex: ['a','b'] -> [True, True, False]
    """
    phase_list = [False, False, False]
    for p in phase_char_lst:
        phase_list[get_phase_idx(p)] = True
    return phase_list

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

def get_Z(dss_data: Any, phase_list : List[bool], fz_mult: float ) -> Iterable:
    """
    helper function to get the Z matrix from dss.lines.to_dataframe()
    Returns an ndarray.
    """
    num_phases = phase_list.count(True)
    if num_phases == 0:
        print(dss_data)
    RM = np.asarray(dss_data['RMatrix'])
    XM = np.asarray(dss_data['XMatrix'])
    ZM = fz_mult * (RM + 1j*XM)
    ZM = np.reshape(ZM, (ZM.shape[0]//num_phases, num_phases))  # reshape
    # pad the Z matrix
    return pad_phases(ZM, (3,3), phase_list)

def pad_phases(matrix: np.ndarray, shape: tuple, phases: List[bool]) -> Iterable:
    """
    Helper method to reshape input matrix and set values set to 0
    for phases set to FALSE in phases tuple.
    Input:
        matrix: an nd array
        shape: a 2-element tuple indicating the output matrix's shape
        phases: a tuple of booleans corresponding to phases to set on this matrix (A: T/F, B: T/F, C: T/F)
    Output:
        matrix reshaped with original values, and
        with 0's for all row/column indices corresponding to phases set to FALSE
    """
    # make the return matrix matrix
    ret_mat = np.zeros(shape, dtype=complex)
    vals = iter(matrix.flatten())
    for out_idx in range(shape[0]):
        if len(shape) == 2:
            for col_idx in range(shape[1]):
                if phases[out_idx] and phases[col_idx]:
                    try:
                        ret_mat[out_idx][col_idx] = next(vals)
                    except StopIteration:
                        ("Cannot pad matrix.")
        else:
            try:
                ret_mat[out_idx] = next(vals)
            except StopIteration:
                ("Cannot pad matrix.")
    return ret_mat

def mask_phases(matrix: Iterable, shape: tuple, phases: List[bool]) -> Iterable:
    """
    Zeroes out values in input matrix for phases set to FALSE in the phases tuple.
    Input:
        matrix: a 3x3 ndarray
        phases: a tuple of booleans corresponding to phases to set on this matrix (A: T/F, B: T/F, C: T/F)
    Output:
        input matrix with 0's for all row/column indices corresponding to phases set to FALSE
    """
    phase_matrix = np.zeros(shape, dtype=complex)
    for out_idx in range(shape[0]):
        if len(shape) == 2:
            for col_idx in range(shape[1]):
                if phases[out_idx] and phases[col_idx]:
                    phase_matrix[out_idx][col_idx] = 1
        else:
            if phases[out_idx] :
                phase_matrix[out_idx] = 1
    masked = np.multiply(matrix, phase_matrix)
    # change all NaN's to 0
    return np.nan_to_num(masked)
