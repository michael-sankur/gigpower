# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Object model for networks
from typing import Any, List, Dict, Iterable, Tuple
import numpy as np
import opendssdirect as dss
import sys
import pandas as pd

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
        self.loads = dict()
        self.capacitors = dict()
        self.adj = dict()
        self.Vbase = None
        self.Sbase = None
        self.Ibase = None
        self.Zbase = None
        self.num_phases = num_phases
        if dss_fp:
            init_from_dss(self, dss_fp)

    def get_nodes(self) -> List[Any] :
        """return a list of network Node objects"""
        return self.nodes.values()

    def get_lines(self) -> List[Any] :
        """return a list of network Line objects"""
        return self.lines.values()

    def __str__(self) -> str:
        """return something informative"""
        df_dict = self.to_data_frames()
        return_str = ''
        for k,v in df_dict.items():
            if k == 'Adj List':
                return_str += '\nAdjacency List:\n'
                for n,l in v.items():
                    return_str += f"{n}: {', '.join(l)}\n"
            else:
                return_str += f"\n\n{k}\n{v}\n"
        return return_str

    def to_dataframes(self) -> Dict[str,Any]:
        """
        returns a dictionary of data frames
        Except for adjacency matrix, which is returned as a list of lists
        """
        # TODO: could define a network object superclass that creates a dataframe from dict
        base_df = pd.DataFrame([self.Vbase, self.Sbase, self.Ibase, self.Zbase, self.num_phases], ['Vbase', 'Sbase', 'Ibase', 'Zbase', 'num_phases']).transpose()
        nodes_df = pd.DataFrame.from_dict( {node.name: node.to_series() for node in self.nodes.values()}).transpose()
        lines_df = pd.DataFrame.from_dict( { str(line.key):line.to_series()  for line in self.lines.values()}).transpose()
        loads_df = pd.DataFrame.from_dict( { load.name:load.to_series()  for load in self.loads.values()}).transpose()
        return({'Base': base_df, 'Nodes': nodes_df, 'Lines': lines_df,'Loads': loads_df, 'Adj List': self.adj})


class Node:
    series_index = ['name', 'phases', 'load', 'controller']
    def __init__(self, name: str = ''):
        self.name = name
        self.phases = (False, False, False)
        self.parent = None # pointer to parent Node. only one parent for radial networks
        self.loads = [] # list of Nodes
        self.caps = [] # list of Capacitors
        self.controller = None
    def __str__(self):
        return f"{self.name}, {self.phases}"
    def to_series(self):
        data = [self.name, self.phases, self.loads, self.controller]
        return pd.Series(data, self.series_index)

class Line:
    series_index = ['(tx,rx)', 'name', 'phases', 'config', 'length', 'FZpu']
    # TODO: might be helpful to include a list of pointers to all Lines in the class, and do the same for Node, etc.
    # see: http://effbot.org/pyfaq/how-do-i-get-a-list-of-all-instances-of-a-given-class.htm
    def __init__(self, key: Tuple[str] = None, name: str = ''):
        self.key = key # tuple of (txnode_name, rxnode_name)
        self.name = name # string, the name DSS uses to refer to this line
        self.phases = (False, False, False)
        self.config = None
        self.length = None
        self.FZpu = np.zeros((3,3), dtype = 'complex')
    def __str__(self):
        return str(self.key)
    def to_series(self):
        data = [self.key, self.name, self.phases, self.config, self.length, self.FZpu]
        return pd.Series(data, self.series_index)



class Load:
    series_index = ['name', 'conn', 'phases', 'type', 'aPQ', 'aI', 'ppu', 'qpu', 'spu']
    def __init__(self, name: str = ''):
        self.name = name
        self.conn = None
        self.phases = (False, False, False)
        self.type = None
        self.aPQ = None
        self.aI = None
        self.ppu = np.zeros((3, 3), dtype='complex')
        self.qpu = np.zeros((3, 3), dtype='complex')
        self.spu = np.zeros((3, 3), dtype='complex')
    def __str__(self):
        return self.name
    def to_series(self):
        data = [self.name, self.conn, self.phases, self.type, self.aPQ, self.aI, self.ppu, self.qpu, self.spu]
        return pd.Series(data, self.series_index)

class Controller:
    def __init__(self, name: str = ''):
        self.node = None
        self.name = None
        self.phases = (False, False, False)
        self.wpu = np.zeroes(3)
        self.wmaxpu = np.zeroes(3)
        self.fes = np.zeroes(3)
        self.hpfes = np.zeroes(3)
        self.lpfes = np.zeroes(3)
        self.kintes = np.zeroes(3)
    def __str__(self):
        return self.name

class Capacitor:
    def __init__(self, name: str = ''):
        self.name = name
        self.phases = (False, False, False)
        self.conn = None
        self.cappu = np.zeros((3, 3), dtype='float')
    def __str__(self):
        return self.name


