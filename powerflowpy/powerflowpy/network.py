# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Object model for networks
from typing import Any, List, Dict, Iterable, Tuple
import numpy as np
import opendssdirect as dss
import sys
import scipy as sp
from utils import init_from_dss

#TODO: determine appropriate precision for attributes. Right now everything is either floats or ints.
# Decimal has alterable precision (defaults to 28 places), and numpy.float128 can have 64-bit precision
# See: https://stackoverflow.com/questions/6663272/double-precision-floating-values-in-python

class Network:
    def __init__(self, dss_fp: str = ''):
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

class Node:
    def __init__(self, name: str = '', num_phases:int = 3, ):
        self.name = name
        self.phases = [0] * num_phases
        self.load = None
        self.controller = None

class Line:
    def __init__(self, key: tuple[str] = None, num_phases:int = 3):
        self.key = key
        self.phases = [0] * num_phases
        self.config = None
        self.length = None
        self.FZpu = np.zeroes((3,3), dtype = 'float')
        self.FRpu = np.zeroes((3, 3), dtype='float')
        self.FXpu = np.zeroes((3, 3), dtype='float')
        self.FYpu = np.zeroes((3, 3), dtype='float')
        self.FGpu = np.zeroes((3, 3), dtype='float')
        self.FBpu = np.zeroes((3, 3), dtype='float')
    
class Load:
    def __init__(self, name: str = '', num_phases=3):
        self.name = name
        self.conn = None
        self.phases = [0] * num_phases
        self.type = None
        self.aPQ = None
        self.aI = None
        self.ppu = None
        self.spu = None

class Controller:
    def __init__(self, name: str = '', num_phases:int = 3):
        self.node = None
        self.phases = [0] * num_phases
        self.wmaxpu = [0] * num_phases
        self.fes = [0] * num_phases
        self.hpfes = [0] * num_phases
        self.lpfes = [0] * num_phases
        self.kintes = [0] * num_phases

class Capacitor:
    def __init__(self, name: str = '', num_phases:int = 3):
        self.name = name
        self.phases = [0] * num_phases
        self.conn = None
        self.cappu = np.zeroes((3, 3), dtype='float')

