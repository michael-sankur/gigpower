# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Implement FBS to solve Power Flow for a radial distribution network.
import sys
from typing import Iterable, List
from utils import is_acyclic, init_from_dss
from network import *

def fbs(dss_fp) -> None:
    if not is_acyclic(dss_fp):
        raise ValueError('Not a radial network.')
    network = init_from_dss(dss_fp)
    # TODO: write fpb!

def _update(network: Network, source_node: Node, target_node: Node) -> None:
    pass

def _topo_sort(network: Network) -> List:
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
    topo_order = [None] * clock  # array to store topo order
    nodes = { node[name]: 'new' for node in network.nodes } # initialize all nodes' statuses to 'New'

    # top level call to _topo_sort_dfs for each node
    #TODO: also return a list of connected commponents, in case we need to know
    #about isolated nodes or the network is a forest of trees
    for node_name,status in nodes.items():
        if status is 'new':
            dfs_stack = [node_name] # use python's list as a stack
            while dfs_stack:
                curr = dfs_stack.pop()
                nodes[curr] = 'active'
                clock -= 1 # decrement clock for a reverse post-order over the vertices
                topo_order[clock] = curr
                for child in network.adj[curr]:
                    if status[child] is 'new':
                        dfs_stack.append(child)
                    elif status[child] is 'active':
                    # child's status is 'active', meaning we have found a back edge
                    # and there is a cycle from child to curr
                        raise ValueError('Network contains a cycle.')
                    # implicit else: child's status is 'finished', so we already assigned the child to topo_order
                status[curr] = 'finished'
            # when the stack is exhausted, we have finished a traversal.

    return topo_order







