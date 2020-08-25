# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Implement FBS to solve Power Flow for a radial distribution network.
import sys
from typing import Iterable
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
    Adopted from Algorithms by Jeff Erickson, Ch. 6.3
    """

    for node in Network._nodes: # use the shadow set of nodes to store topo sort info



