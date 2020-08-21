# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Object model for networks
from typing import Any, List, Dict, Iterable
import numpy as np
import opendssdirect as dss
import sys
import scipy as sp
from utils import init_from_dss

class Network:
    def __init__(self, dss_fp: str = None):
        """ 
        Initialize a Network instance. 
        Optional dss_file is the full path to a dss file.
        """
        self.nodes = dict()
        self.lines = dict()
        self.adj = dict()
        self.Vbase = None
        self.Sbase = None
        self.Ibase = None
        self.Zbase = None
        self.units = None
        self._nodes = dict() # stores topology info from dfs
        if dss_fp:
            init_from_dss(self, dss_fp)

    def topo_sort(self) -> List:
        """returns list of Nodes in topological sort order""""
        #TODO: handle cases of multiple valid sorts
        if self.has_cycle():
            raise ValueError('Network contains a cycle.')
        pass

    def has_cycle(self) -> bool:
        """returns true if Network has a cycle"""
        pass

    def _dfs(self)-> None:
        """traverse network by dfs, store traversal information on self._nodes"""
        pass
    
    
    
