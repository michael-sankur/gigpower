# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Implement FBS to solve Power Flow for a radial distribution network.
import sys
import time
import csv
from typing import Iterable, List, Dict
from . utils import init_from_dss, mask_phases
from . network import *
from . solution import *

def fbs(dss_file) -> List[float]:

    # MAP NETWORK FROM DSS FILE-------------------------------------------------
    t1 = time.perf_counter()
    network = init_from_dss(dss_file)
    t2 = time.perf_counter()
    init_from_dss_time = (t2 - t1) * 1000

    # TOPO SORT-----------------------------------------------------------------
    t3 = time.perf_counter()
    topo_order = topo_sort(network)
    t4 = time.perf_counter()
    topo_sort_time = (t4 - t3) * 1000

    # INITIALIZE SOLUTION OBJECT------------------------------------------------
    t5 = time.perf_counter()
    solution = Solution(network)
    solution.root = network.nodes.get(topo_order[0])  # keep a pointer to the root node
    solution.tolerance = abs( (solution.Vref[1]) * 10**-9) # set tolerance with phase B reference voltage
    converged = max(abs(solution.Vtest - solution.Vref)) <= solution.tolerance
    # init V.root to Vref
    solution.V[solution.root.name] = np.copy(solution.Vref)
    t6 = time.perf_counter()
    solution_init = (t6 - t5) * 1000

    while not converged:
        fwd_sweep_total = 0
        pre_bwd_sweep_load_update_total = 0
        bwd_sweep_total = 0
        convergence_check_total = 0

    # FORWARD SWEEP ------------------------------------------------------------
        t7 = time.perf_counter()
        # FORWARD SWEEP: for node in topo_order:
        for node_name in topo_order:
            node = network.nodes.get(node_name)
            children = network.adj[node_name]
            for child_name in children:
                child = network.nodes.get(child_name)
                line_out = network.lines[(node_name, child_name)]
                solution.update_voltage_forward(network, node, child)
        t8 = time.perf_counter()
        fwd_sweep_total += (t8-t7) * 1000

    # PRE- BACKWARD SWEEP, UPDATE VOLTAGE DEPENDENT LOAD AT ALL NODES-----------
        for node in network.get_nodes():
            solution.update_voltage_dependent_load(node) # update s at all nodes
        t9 = time.perf_counter()
        pre_bwd_sweep_load_update_total += (t9 - t8) * 1000

    # BACKWARD SWEEP -----------------------------------------------------------
        # BACKWARD SWEEP: for node in reverse topo_order:
        for node_name in reversed(topo_order):
            node = network.nodes.get(node_name)
            if node.parent: # if this is a terminal node or junction node (not the root)
                solution.update_voltage_dependent_load(
                    node)  # update s at this node
                solution.update_voltage_dependent_load(
                    node.parent)  # update s at node's parent
                line_in = network.lines.get((node.parent.name, node.name))
                solution.update_current(network, line_in) # update current segment
                solution.update_voltage_backward(
                    network, node)  # update voltage at parent
        t10 = time.perf_counter()
        bwd_sweep_total += (t10 - t9) * 1000

    # CHECK CONVERGENCE---------------------------------------------------------
        solution.iterations += 1
        # set Vtest to the root's voltage
        solution.Vtest = solution.V[solution.root.name]
        solution.diff = max(abs(solution.Vtest - solution.Vref))

        #check convergence
        converged = max(abs(solution.Vtest - solution.Vref)
                        ) <= solution.tolerance
        # if we're looping again, set V.root to V.ref
        if not converged:
            solution.V[solution.root.name] = np.copy(solution.Vref)
        t11 = time.perf_counter()

        convergence_check_total += (t11 - t10) * 1000

    # FINAL CALCULATIONS -------------------------------------------------------
    solution.calc_sV()
    solution.calc_S()
    solution.calc_Inode()
    t12 = time.perf_counter()
    final_calcs = (t12 - t11) * 1000

    return [ init_from_dss_time, topo_sort_time, solution_init, fwd_sweep_total, pre_bwd_sweep_load_update_total, bwd_sweep_total, convergence_check_total, final_calcs ]


def topo_sort(network: Network) -> List:
    """
    Returns a list of Nodes in a topological order.
    Raises error if it finds a cycle
    Adapted from 'Algorithms' by Jeff Erickson, 1st ed 2019, Ch. 6.3, Figure 6.9,
    'Explicit topological sort': http://jeffe.cs.illinois.edu/teaching/algorithms/book/06-dfs.pdf
    """
    # store each node's status. 'New' means it has not been traversed,
    # 'Active' means it is being explored.
    clock = len(network.nodes) # countdown from total number of nodes.
    # clock will be used as an index to topo_order. the strategy is to insert
    # nodes in the topo_order array in reverse order of dfs finishing times.
    topo_order = [None] * clock # array to store topo order
    nodes = { node_name: 'new' for node_name in network.nodes.keys()} # initialize all nodes' statuses to 'New'

    # top level call to _topo_sort_dfs for each node
    for node_name,status in nodes.items():
        if status is 'new':
            clock = topo_sort_dfs(node_name, nodes, network.adj, clock, topo_order)
    return topo_order

def topo_sort_dfs(start_node:str, node_status: Dict[str,str], adj_matrix: Iterable, clock: int, topo_order: List) -> int:
    node_status[start_node] = 'active'
    for child in adj_matrix[start_node]:
        if node_status[child] is 'new':
            clock = topo_sort_dfs(child, node_status, adj_matrix, clock, topo_order)
        elif node_status[child] is 'active':
            raise ValueError('Network contains a cycle.')
    node_status[start_node] = 'finished'
    topo_order[clock-1] = start_node
    return clock - 1


