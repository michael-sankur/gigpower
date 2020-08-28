# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Implement FBS to solve Power Flow for a radial distribution network.
import sys
from typing import Iterable, List, Dict
from . utils import init_from_dss
from . network import *

def fbs(dss_fp) -> None:
    network = init_from_dss(dss_fp)
    topo_order = topo_sort(network)
    solution = Solution(network)
    #TODO: figure out how best to set tolerance
    solution.tolerance = abs( (solution.Vref[1]) * 10**-9) # set tolerance at phase B

    converged = max(abs(Vtest - Vref)) >= solution.tol
    while not converged:
        # for node in topo_order:
            # forward sweep
        converged = max(abs(Vtest - Vref)) >= solution.tol
        # for node in reverse topo_order:
            # backward sweep
    return Solution

def update(network: Network, source_node: Node, target_node: Node) -> None:
    pass

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


class Solution:

    def __init__(self, network: Network, tol: float = -1, max_iter: int = -1) -> None:
        self.iterations = -1  # stores number of iterations of FBS until convergence
        self.Vtest = np.zeros(3, dtype='complex')
        self.Vref = network.Vbase * np.ones(3, dtype='complex')
        self.solved_nodes = dict()
        self.solved_lines = dict()
        self.tolerance = -1 # stores the tolerance
        self.diff = -1  # stores the final value of Vtest - Vref at convergence

        """ Set up voltages and currents for all nodes """
        for node in network.get_nodes():
            node_dict = dict()
            self.solved_nodes[node.name] = node_dict
            # initialize a zeroed out array for each solution var
            node_dict['V'] = network.Vbase * np.ones(3, dtype='complex')
            node_dict['Inode'] = np.zeros(3, dtype='complex')
            node_dict['S'] =
        for line in network.get_lines():
            line_dict = dict()
            self.solved_lines[line.name] = line_dict
            line_dict['I'] = np.zeros(3s, dtype='complex')

    def __str__(self):
        return '\n'.join(
            [
                f'itertations to convergence: {self.iterations}',
                f'tolerance at convergence: {self.tolerance}',
                f'solution: \n {self.solved_nodes}'
            ]
        )





