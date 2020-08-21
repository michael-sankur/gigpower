# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Implement FBS to solve Power Flow for a radial distribution network. 
import sys
from typing import Iterable
from utils import is_acyclic, init_from_dss

def fbs(dss_fp) -> None:
    if not is_acyclic(dss_fp):
        raise ValueError('Not a radial network.')
    network = init_from_dss(dss_fp)
    # TODO: write fpb!

def _update(network: Network, source_node: Node, target_node: Node) -> None:
    pass

def _topo_sort(network: Network) -> List:
    """returns list of Nodes in topological sort order"""
    #TODO: handle cases of multiple valid sorts
    if self.has_cycle():
        raise ValueError('Network contains a cycle.')
    pass

def has_cycle(network:Network) -> bool:
    """returns true if Network has a cycle"""
    pass

def _dfs(network:Network ) -> None:
    """traverse network by dfs, store traversal information on self._nodes"""
    pass

