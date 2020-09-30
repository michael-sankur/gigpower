# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Implement FBS to solve Power Flow for a radial distribution network.
import sys
from typing import Iterable, List, Dict
from . utils import init_from_dss, mask_phases
from . network import *
from . solution import *

def fbs(network, dss_fp = True) -> None:
    if dss_fp:
        network = init_from_dss(network)
    solution = Solution(network)
    topo_order = topo_sort(network)
    solution.root = network.nodes.get(topo_order[0])  # keep a pointer to the root node
    #TODO: Make a better 'solution.set_tolerance(ref_node, error)' method
    solution.tolerance = abs( (solution.Vref[1]) * 10**-9) # set tolerance with phase B reference voltage
    converged = max(abs(solution.Vtest - solution.Vref)) <= solution.tolerance
    while not converged:
        # set V.root to Vref
        solution.V[solution.root.name] = np.copy(solution.Vref)

        # FORWARD SWEEP: for node in topo_order:
        for node_name in topo_order:
            node = network.nodes.get(node_name)
            children = network.adj[node_name]
            for child_name in children:
                child = network.nodes.get(child_name)
                line_out = network.lines[(node_name, child_name)]
                solution.update_voltage_forward(network, node, child)

        for node in network.get_nodes():
            solution.update_voltage_dependent_load(node) # update s at all nodes

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

        solution.iterations += 1
        # set Vtest to the root's voltage
        solution.Vtest = solution.V[solution.root.name]
        solution.diff = max(abs(solution.Vtest - solution.Vref))
        #check convergence
        converged = max(abs(solution.Vtest - solution.Vref)
                        ) <= solution.tolerance

    # # final calculations
    # solution.calc_sV()
    # solution.calc_S()
    # # TODO: check that node_sV = node_Srx - sum(all line.Stx for all node.outgoing_lines)
    # solution.calc_Inode()
    return solution

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
    #TODO: also return a list of connected commponents, in case we need to know
    #about isolated nodes or the network is a forest of trees
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


