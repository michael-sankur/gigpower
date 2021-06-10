from typing import List, Union, Dict
import numpy as np
import pandas as pd


def parse_dss_bus_name(dss_bus_name: str, sep='.') -> str:
    """
    Given a bus name string that may include phase from opendss, returns just
    the busname.
    Assumes that dss appends bus names with phases, separated by '.'
    Ex: 'sourcebus.1.2.3' -> 'sourcebus'
    """
    return dss_bus_name.split(sep)[0]


def parse_dss_phases(dss_str: str, sep='.') -> List[str]:
    """
    Given a bus name string that may include phase from opendss, returns a list
    of chars representing phases.
    Assumes that dss appends bus names with phases, separated by '.'
    Ex: 'sourcebus.2.3' -> ['2', '3']
    If there are no phase numbers, assumes that all phases are present on the bus.
    """
    if sep in dss_str:
        return [str(ph) for ph in dss_str.split(sep)[1:]]
    else:
        return ['1', '2', '3']


def parse_phases(phase_list: List[Union[str, int]]) -> np.ndarray:
    """
    helper function to return a list of phases represented as ints or
    strings into a list of booleans
    ex: ['1', '3'] -> [True, False True]
    ex: ['a','b'] -> [True, True, False]
    """
    phase_bools = [False, False, False]
    for p in phase_list:
        phase_bools[get_phase_idx(p)] = True
    return np.asarray(phase_bools)


def parse_phase_matrix(phase_char_lst: List[str]) -> np.ndarray:
    """
    n x 3 phase matrix of 1's where phases are present, 0's otherwise
    """
    return np.asarray([1 if ph else 0 for ph in parse_phases(phase_char_lst)])


def get_phase_idx(ph: Union[str, int]) -> int:
    """
    helper function to turn a phase letter into an index, where 'a' = 0
    """
    if ph in ['a', 'b', 'c']:
        return ord(ph.lower()) - ord('a')
    elif ph in ['1', '2', '3']:
        return int(ph) - 1
    elif ph in range(1, 4):
        return ph - 1
    else:
        raise ValueError(f'Invalid argument for get_phase_idx {ph}')


def set_zip_values_dss(dss, zipv: List):
    """sets zip values on a dss object by setting the dss.Loads.ZipV"""
    # array mapping: [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min votlage pu]
    for load_name in dss.Loads.AllNames():
        dss.Loads.Name(load_name)
        dss.Loads.Model(8)
        dss.Loads.ZipV(zipv)
        dss.Loads.ZipV()


def pad_phases(matrix: np.ndarray, shape: tuple, phases: np.ndarray) -> np.ndarray:
    """
    Helper method to reshape input matrix according to a phases array
    Input:
        matrix: an nd array
        shape: a 2-element tuple indicating the output matrix's shape
        phases: 3 x 1 nd array with ones indicating phases present, 0's otherwshie
    Output:
        matrix reshaped with original values, and
        with 0's for all row/column indices corresponding to non-existent phases
    """
    # make the return matrix matrix
    ret_mat = np.zeros(shape, dtype=matrix.dtype)
    vals = iter(matrix.flatten())
    for out_idx in range(shape[0]):
        if len(shape) == 2:
            for col_idx in range(shape[1]):
                if phases[out_idx] == 1 and phases[col_idx] == 1:
                    try:
                        ret_mat[out_idx][col_idx] = next(vals)
                    except StopIteration:
                        ("Cannot pad matrix.")
        else:
            try:
                if phases[out_idx] == 1:
                    ret_mat[out_idx] = next(vals)
            except StopIteration:
                ("Cannot pad matrix.")
    return ret_mat


def mask_phases(matrix: np.ndarray, shape: tuple, phases: np.ndarray) -> np.ndarray:
    """
    Zeroes out values in input matrix for phases set to ZERO in the phase matrix
    Input:
        matrix: a 3x3 ndarray
        phases: nd array with 1's where phases exist, zeros otherwise
    Output:
        input matrix with 0's for all row/column indices corresponding 
        to phases that are 0 in the phase matrix
    """
    out_matrix = np.zeros(shape, dtype=matrix.dtype)
    for out_idx in range(shape[0]):
        if len(shape) == 2:
            for col_idx in range(shape[1]):
                if phases[out_idx] == 1 and phases[col_idx] == 1:
                    out_matrix[out_idx][col_idx] = 1
        else:
            if phases[out_idx]:
                out_matrix[out_idx] = 1
    masked = np.multiply(matrix, out_matrix)
    # change all NaN's to 0
    return np.nan_to_num(masked)



def topo_sort(bus_names: List[str], adj_matrix: Dict) -> List:
    """
    Returns a list of Busesin a topological order.
    Raises error if it finds a cycle
    Adapted from 'Algorithms' by Jeff Erickson, 1st ed 2019, Ch. 6.3, Figure 6.9,
    'Explicit topological sort': http://jeffe.cs.illinois.edu/teaching/algorithms/book/06-dfs.pdf
    """
    # store each node's status. 'New' means it has not been traversed,
    # 'Active' means it is being explored.
    clock = len(bus_names)  # countdown from total number of nodes.
    # clock will be used as an index to topo_order. the strategy is to insert
    # buses in the topo_order array in reverse order of dfs finishing times.
    topo_order = [None] * clock  # array to store topo order
    # initialize all nodes' statuses to 'New'
    buses = {bus_name: 'new' for bus_name in bus_names}

    # top level call to _topo_sort_dfs for each node
    # TODO: also return a list of connected commponents, in case we need to know
    # about isolated nodes or the network is a forest of trees
    for bus_name, status in buses.items():
        if status == 'new':
            clock = topo_sort_dfs(
                bus_name, buses, adj_matrix, clock, topo_order)
    return topo_order

def topo_sort_dfs(start_node: str, node_status: Dict[str, str],
                    adj_matrix: Dict[str, List[str]], clock: int, topo_order: List) -> int:
    node_status[start_node] = 'active'
    if start_node in adj_matrix: 
        for child in adj_matrix[start_node]:
            if node_status[child] == 'new':
                clock = topo_sort_dfs(child, node_status,
                                            adj_matrix, clock, topo_order)
            elif node_status[child] == 'active':
                raise ValueError('Not radial. Contains a cycle.')
    node_status[start_node] = 'finished'
    topo_order[clock-1] = start_node
    return clock - 1


def get_reverse_adj(adj: Dict[str, List]):
    """ Returns the reverse adjacency matrix"""
    reverse = {}
    for node, neighbors in adj.items():
        if node not in reverse:
            reverse[node] = []
        for v in neighbors:
            try:
                reverse[v].append(node)
            except KeyError:
                reverse[v] = [node]
    return reverse


def get_nominal_bus_powers(dss) -> pd.DataFrame:
    """
    Interface directly with dss object to get nominal bus powers of current
    circuit. Used for testing.
    """
    bus_names = dss.Circuit.AllBusNames()
    bus_idx_dict = {b: i for i, b in enumerate(bus_names)}
    data = np.zeros((len(bus_names), 3), dtype=complex)

    # add all load powers to load data, based on bus index
    for load in dss.Loads.AllNames():
        dss.Loads.Name(load)
        bus_name = dss.CktElement.BusNames()[0]
        bus_name, bus_phase = bus_name.split('.')[0], bus_name.split('.')[1:]
        if len(bus_phase) == 0:
            bus_phase.extend(['1', '2', '3'])
        bus_idx = bus_idx_dict[bus_name]
        phases = parse_phases(list(bus_phase))
        dist_real = dss.Loads.kW() / len(bus_phase)
        dist_imag = dss.Loads.kvar() / len(bus_phase)
        data[bus_idx, phases] += (dist_real + 1j * dist_imag)

    # add all capacitor powers to load data, based on bus index
    for cap in dss.Capacitors.AllNames():
        dss.Capacitors.Name(cap)
        bus_name = dss.CktElement.BusNames()[0]
        bus_name, bus_phase = bus_name.split('.')[0], bus_name.split('.')[1:]
        if len(bus_phase) == 0:
            bus_phase.extend(['1', '2', '3'])
        bus_idx = bus_idx_dict[bus_name]
        phases = parse_phases(list(bus_phase))
        real, imag = 0, dss.Capacitors.kvar() / (phases == True).sum()
        data[bus_idx, phases] += (real - 1j * imag)

    return pd.DataFrame(data, bus_names, ['A', 'B', 'C'])
