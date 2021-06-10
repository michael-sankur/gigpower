import numpy as np # type: ignore
import opendssdirect as dss # type: ignore
from typing import Iterable, List, Any
from . network import Network, Node, Line, Load, Capacitor, Transformer, VoltageRegulator

ZIPV = [0.10, 0.05, 0.85, 0.10, 0.05, 0.85, 0.80]
# ZIPV = [0, 0, 1, 0, 0, 1, 0.75]
# ZIPV = [1, 0, 0, 1, 0, 0, 0.75] # constance impedance
ZIPV_caps = [1, 0, 0, 1, 0, 0]  # for reference, unused


def init_from_dss(dss_fp: str) -> Network:
    """define a Network attributes from a dss file"""
    network = Network()
    dss.run_command('Redirect ' + dss_fp)
    # solve, get Base values
    dss.Solution.Solve()
    # set the load model and zip values
    set_zip_values(dss, ZIPV)
    # solve again
    dss.Solution.Solve()
    get_nodes_from_dss(network, dss)
    get_lines_from_dss(network, dss)
    get_loads_from_dss(network, dss)
    get_caps_from_dss(network, dss)

    # note: map voltage regulators before mapping transformers
    # so that we only map transformers that are not voltage regulators
    get_voltage_regulators_from_dss(network, dss)
    get_transformers_from_dss(network, dss)

    return network


def set_zip_values(dss:Any, zipv: List) -> None:
    """sets custom zip values in dss by setting the dss.Loads.zipv() array."""
    # array mapping: [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min votlage pu]
    for load_name in dss.Loads.AllNames():
        dss.Loads.Name(load_name)
        dss.Loads.Model(8)
        dss.Loads.ZipV(zipv)
        dss.Loads.ZipV()


def get_nodes_from_dss(network: Network, dss: Any) -> None:
    # make Nodes
    for node_name in dss.Circuit.AllNodeNames():
        name, phase = node_name.split('.')
        # get the node corresponding to this name, or make a new one
        if name not in network.nodes.keys():
            # initialize node
            node = Node(name)
            node.phases = []  # replace default tuple with empty list so we can mutate it
            dss.Circuit.SetActiveBus(name)
            node.Vbase = dss.Bus.kVBase() * 1000
            node.Ibase = node.Sbase/node.Vbase
            node.Zbase = node.Vbase/node.Ibase
            network.nodes[name] = node  # add node to network
            # Add node to adjacency list. note: this means that every node has an entry.
            # Nodes with no children will have an empty list.
            network.adj[name] = []
        node = network.nodes[name]

        # store phases as characters in phase list for now
        node.phases.append(phase)

    # iterate through nodes to parse phase lists
    for node in network.get_nodes():
        node.phases = parse_phases(node.phases)  # type: ignore

    # save a bus name: index dict on the network
    network.bus_idx_dict = {b: i for i, b in enumerate(dss.Circuit.AllBusNames())}


def get_lines_from_dss(network: Network, dss: Any) -> None:
    # make Lines
    all_lines_data = dss.utils.lines_to_dataframe().transpose()  # get dss line data indexed by line_code
    line_codes = all_lines_data.keys()
    # TODO: set line base values based on node Vbase
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

        line = Line(network, (tx, rx), line_code)  # initialize line.
        # Line constructor will add the line to the adjacency list, and set the rx.parent to tx
        # set phases according to tx
        line.phases = parse_phases(tx_phases)
        tx_node = network.nodes.get(tx)
        # parse line attributes from dss line data
        line.name = line_code
        line.length = line_data['Length']
        fz_mult = 1 / tx_node.Zbase * line.length  # use tx node's Z base
        line.FZpu = get_Z(line_data, line.phases, fz_mult)


def get_loads_from_dss(network: Network, dss: Any) -> None:
    # make Loads
    all_loads_data = dss.utils.loads_to_dataframe().transpose()
    load_names = all_loads_data.keys()
    for load_name in load_names:
        load_data = all_loads_data[load_name]
        dss.Loads.Name(load_name)
        bus_phase = dss.CktElement.BusNames()[0].split('.')
        if len(bus_phase) == 1:
            bus_phase.extend(['1', '2', '3'])
        node_name, phase_chars = bus_phase[0], bus_phase[1:]
        try:
            node = network.nodes[node_name]
        except KeyError:
            print(f"This load's node has not been defined. Load name: {load_name}, Node name: {node_name}")
        load = Load(load_name)
        load.phases = parse_phases(list(phase_chars))
        # save kw and kvar
        load.kw = load_data['kW']
        load.kvar = load_data['kvar']
        load.conn = 'delta' if load_data['IsDelta'] else 'wye'

        # divide ppu and qpu by number of phases
        load.ppu = load.kw / 1000  / load.phases.count(True)
        load.qpu = load.kvar / 1000  / load.phases.count(True)
        # set aPQ, aI, aZ
        if load_data['ZipV']:
            load.zipV = load_data['ZipV']
            # array mapping: [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min votlage pu]
            load.aZ_p = load.zipV[0]
            load.aI_p = load.zipV[1]
            load.aPQ_p = load.zipV[2]
            load.aZ_q = load.zipV[3]
            load.aI_q = load.zipV[4]
            load.aPQ_q = load.zipV[5]
        else:
            load.aPQ_p, load.aI_p, load.aZ_p = 1, 0, 0
            load.aPQ_q, load.aI_q, load.aZ_q = 1, 0, 0
        load.spu = load.ppu + 1j*load.qpu

        # add a pointer to this load to the network
        network.loads[load_name] = load
        load.node = node  # assign the node to the load
        node.loads.append(load)  # add this load to its node's load list
        # TODO: set load.type


def get_caps_from_dss(network: Network, dss: Any) -> None:
    # make Capacitors
    all_cap_data = dss.utils.capacitors_to_dataframe().transpose()
    cap_names = all_cap_data.keys()
    for cap_name in cap_names:
        if cap_name != '':
            cap_data = all_cap_data[cap_name]
            dss.Capacitors.Name(cap_name)
            bus_phase = dss.CktElement.BusNames()[0].split('.')
            if len(bus_phase) == 1:
                bus_phase.extend(['1', '2', '3'])
            node_name, phase_chars = bus_phase[0], bus_phase[1:]
            try:
                node = network.nodes[node_name]
            except KeyError:
                print(
                    f"This cap's node has not been defined. Cap name: {cap_name}, Node name: {node_name}")
            cap = Capacitor(cap_name)
            cap.phases = parse_phases(list(phase_chars))
            cap.conn = 'delta' if cap_data['IsDelta'] else 'wye'
            cappu = cap_data['kvar'] * 1000 / network.Sbase / cap.phases.count(True)
            cap.cappu = np.asarray([cappu if phase else 0 for phase in cap.phases])
            network.capacitors[cap_name] = cap  # add a pointer to this cap to the network
            node.capacitors.append(cap)  # add capacitor to it's node's cap list
    # sum all cappu on each node and store on each node, to avoid re-calculating during
    # fbs.update-voltage-dependent-load()
    for node in network.nodes.values():
        for cap in node.capacitors:
            node.sum_cappu = np.add(node.sum_cappu, cap.cappu)


def get_voltage_regulators_from_dss(network: Network, dss: Any) -> None:
    # get Reg Controls from dss
    regs = dss.RegControls.AllNames()
    for regControl_name in regs:
        dss.RegControls.Name(regControl_name)  # set this reg as active, to get the transformer
        transformer_name = dss.RegControls.Transformer()
        dss.Transformers.Name(transformer_name)  # set this regcontrol's transformer as active
        tx, regControl_node = dss.CktElement.BusNames()  # get upstream, regcontrol nodes
        tx, *phases = tx.split('.')
        regControl_node = regControl_node.split('.')[0]
        voltageReg = VoltageRegulator(network, regControl_name, regControl_node, tx)
        voltageReg.transformer_name = transformer_name
        dss.RegControls.Name(regControl_name)  # set this reg as active, again, to get the tap number
        voltageReg.tap = dss.RegControls.TapNumber()
        voltageReg.get_gamma(dss.RegControls.TapNumber())
        voltageReg.phases = parse_phases(phases)


def get_transformers_from_dss(network: Network, dss: Any) -> None:
    # make Transformers
    all_transformer_data = dss.utils.transformers_to_dataframe().transpose()
    transformers = dss.Transformers.AllNames()
    # exclude the transformers that are paired with voltage regulators
    vr_transformers = [r.transformer_name for r in network.voltageRegulators.values()]
    transformer_names = [n for n in transformers if n not in vr_transformers]

    for transformer_name in transformer_names:
        transformer_data = all_transformer_data[transformer_name]
        dss.Transformers.Name(transformer_name)

        bus1, bus2 = (b.split('.')[0] for b in dss.CktElement.BusNames())
        transformer = Transformer(network, (bus1, bus2), transformer_name, transformer_data['NumWindings'])

        dss.Transformers.Wdg(1)  # upstream side
        if dss.CktElement.NumPhases() == 1: # 1 phase -> LN voltage,  2 or 3 -> LL voltage
            upstream = dss.Transformers.kV()
        else:
            upstream = dss.Transformers.kV() / (3 ** 0.5)

        dss.Transformers.Wdg(2)  # downstream side
        if dss.CktElement.NumPhases() == 1:  # 1 phase -> LN voltage,  2 or 3 -> LL voltage
            downstream = dss.Transformers.kV()
        else:
            downstream = dss.Transformers.kV() / (3 ** 0.5)
        transformer.rated_voltages = (upstream, downstream)

        transformer.phases = network.nodes[bus1].phases
        transformer.conn = 'delta' if transformer_data['IsDelta'] else 'wye'

        transformer.kV = transformer_data['kV']
        transformer.kVA = transformer_data['kVA']


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


def get_Z(dss_data: Any, phase_list: List[bool], fz_mult: float) -> Iterable:
    """
    helper function to get the Z matrix from dss.lines.to_dataframe()
    Returns an ndarray.
    """
    num_phases = phase_list.count(True)
    RM = np.asarray(dss_data['RMatrix'])
    XM = np.asarray(dss_data['XMatrix'])
    ZM = fz_mult * (RM + 1j*XM)
    ZM = np.reshape(ZM, (ZM.shape[0]//num_phases, num_phases))  # reshape
    # pad the Z matrix
    return pad_phases(ZM, (3, 3), phase_list)


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
                if phases[out_idx]:
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
            if phases[out_idx]:
                phase_matrix[out_idx] = 1
    masked = np.multiply(matrix, phase_matrix)
    # change all NaN's to 0
    return np.nan_to_num(masked)


def calc_load_power(node: Node, nodeV: Iterable, zipV: Iterable = None) -> Iterable:
    """
    helper function to calculate load powers over a node, given the solved V value, and zip values
    If zipV is not passed by argument, will default to values saved on the Load objects
    pointed to by the given Node

    node: Node object
    nodeV: 3x1 complex volgaget
    zipV: 6x1 zipV array
    """
    abs_nodeV = abs(nodeV)
    abs_nodeV_sq = np.power(abs(nodeV), 2)
    load_powers = np.zeros(3, dtype=complex)

    for load in node.loads:
        spu_real, spu_imag = load.ppu, load.qpu
        if not zipV:
            aPQ_p, aI_p, aZ_p = load.aPQ_p, load.aI_p, load.aZ_p
            aPQ_q, aI_q, aZ_q = load.aPQ_q, load.aI_q, load.aZ_q
        else:
            aZ_p = zipV[0]
            aI_p = zipV[1]
            aPQ_p = zipV[2]
            aZ_q = zipV[3]
            aI_q = zipV[4]
            aPQ_q = zipV[5]

        for idx, ph in enumerate(load.phases):
            if ph:
                temp1 = aPQ_p + aI_p * abs_nodeV[idx] + aZ_p * abs_nodeV_sq[idx]
                real = temp1 * spu_real

                temp2 = aPQ_q + aI_q * abs_nodeV[idx] + aZ_q * abs_nodeV_sq[idx]
                imag = temp2 * spu_imag

                load_powers[idx] += real + (1j * imag)
    return load_powers


def calc_cap_power(node: Node, nodeV: Iterable) -> Iterable:
    """
    helper function to calculate cap powers over a node, given the solved V value

    node: Node object
    nodeV: 3x1 complex volgage
    """
    abs_nodeV_sq = np.power(abs(nodeV), 2)
    cap_powers = np.zeros(3, dtype=complex)

    for cap in node.capacitors:
        cappu = cap.cappu
        cap.imag = np.zeros(3)
        for idx, ph in enumerate(cap.phases):
            if ph:
                real = 0
                imag = cappu[idx] * abs_nodeV_sq[idx]
                cap.imag[idx] = imag
                cap_powers[idx] += real - (1j * imag)
    return cap_powers


def calc_total_node_power(node: Node, nodeV: Iterable, zipV: Iterable = None) -> Iterable:
    """
    helper function to calculate total node power given a Node, solved V value, and zip values
    If zipV is not passed by argument, will default to values saved on the Load objects
    pointed to by the given Node

    node: Node object
    nodeV: 3x1 complex volgaget
    zipV: 6x1 zipV array
    """
    wpu = np.zeros(3)  # TODO: will be set as argument

    total_powers = calc_load_power(node, nodeV, zipV)
    total_powers += calc_cap_power(node, nodeV)
    total_powers += wpu

    return total_powers