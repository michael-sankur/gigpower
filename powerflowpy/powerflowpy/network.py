# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Object model for networks
from typing import Any, List, Dict, Iterable, Tuple
import numpy as np
import opendssdirect as dss
import sys
import scipy as sp

#TODO: determine appropriate precision for attributes. Right now everything is either floats or ints.
# Decimal has alterable precision (defaults to 28 places), and numpy.float128 can have 64-bit precision
# See: https://stackoverflow.com/questions/6663272/double-precision-floating-values-in-python

class Network:
    def __init__(self, dss_fp: str = '', num_phases = 3):
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
        self.phases = [0] * num_phases
        if dss_fp:
            init_from_dss(self, dss_fp)

    def __str__(self):
        """return something informative"""
        node_str = 'Nodes: ' + ', '.join([k for k  in self.nodes.keys()])
        lines_str = 'Lines: ' + ', '.join([ str(k) for k in self.lines.keys()])
        adj_str = f'Adjacency List: {self.adj}'
        params_str = f'Vbase: {self.Vbase}, Sbase: {self.Sbase}, Ibase: {self.Ibase}, Zbase: {self.Zbase}, phases: {len(self.phases)}, units: {self.units}'
        return '\n'.join([params_str, node_str, lines_str, adj_str]) + '\n'

class Node:
    def __init__(self, name: str = '', num_phases:int = 3, ):
        self.name = name
        self.phases = [0] * num_phases
        self.load = None
        self.controller = None
    def __str__(self):
        return self.name

class Line:
    def __init__(self, key: Tuple[str] = None, num_phases:int = 3):
        self.key = key
        self.phases = [0] * num_phases
        self.config = None
        self.length = None
        self.FZpu = np.zeros((3,3), dtype = 'float')
        self.FRpu = np.zeros((3, 3), dtype='float')
        self.FXpu = np.zeros((3, 3), dtype='float')
        self.FYpu = np.zeros((3, 3), dtype='float')
        self.FGpu = np.zeros((3, 3), dtype='float')
        self.FBpu = np.zeros((3, 3), dtype='float')
    def __str__(self):
        return self.key

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
    def __str__(self):
        return self.name

class Controller:
    def __init__(self, name: str = '', num_phases:int = 3):
        self.node = None
        self.name = None
        self.phases = [0] * num_phases
        self.wmaxpu = [0] * num_phases
        self.fes = [0] * num_phases
        self.hpfes = [0] * num_phases
        self.lpfes = [0] * num_phases
        self.kintes = [0] * num_phases
    def __str__(self):
        return self.name

class Capacitor:
    def __init__(self, name: str = '', num_phases:int = 3):
        self.name = name
        self.phases = [0] * num_phases
        self.conn = None
        self.cappu = np.zeros((3, 3), dtype='float')
    def __str__(self):
        return self.name

class Solution:
    # class variable: list of solution variables
    solution_vars = [ 'V', 'I', 'Inode', 'S', 'sV']

    def __init__(self, nodes: Iterable[Node], iterations:int = -1, tol:float = -1) -> Dict[str, Dict] :
        self.iterations = iterations # number of iterations of FBS until convergence
        self.tolerance = tol # tolerance at convergence
        self.solved_nodes = dict()
        for node in self.nodes:
            # initialize a zeroed out num_phases x 1 array for each solution var
            # for each node
            node_dict = { var : [0] * len(node.phases) for var in Solution.solution_vars }
    def __str__(self):
        return '\n'.join(
            [ 
                f'itertations to convergence: {self.iterations}',
                f'tolerance at convergence: {self.tolerance}',
                f'solution: \n {self.solved_nodes}'        
            ]
        )


