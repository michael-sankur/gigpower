# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 20 August 2020
# Object model for networks
from typing import Any, List, Dict, Tuple, ValuesView, Union
import numpy as np # type: ignore
import pandas as pd # type: ignore

#TODO: determine appropriate precision for attributes. Right now everything is either floats or ints.
# Decimal has alterable precision (defaults to 28 places), and numpy.float128 can have 64-bit precision
# See: https://stackoverflow.com/questions/6663272/double-precision-floating-values-in-python


class Network:
    def __init__(self, num_phases:int = 3) -> None:
        """
        Initialize a Network instance.
        """
        self.nodes: Dict[str,Node] = dict()
        self.lines: Dict[Tuple [str,str], Line] = dict()
        self.loads: Dict[str, Load] = dict()
        self.capacitors: Dict[str, Capacitor] = dict()
        self.transformers: Dict[str, Transformer] = dict()
        self.adj: Dict[str, List] = dict()
        self.Vbase = 0.0
        self.Sbase = 0.0
        self.Ibase = 0.0
        self.Zbase = 0.0
        self.num_phases = num_phases


    def get_nodes(self) -> ValuesView[Any]:
        """return a List-like iterable view of network Node objects"""
        return self.nodes.values()

    def get_lines(self) -> ValuesView[Any]:
        """return a List-like view of network Line objects"""
        return self.lines.values()

    def get_loads(self) -> ValuesView[Any]:
        """return a List-like view network Load objects"""
        return self.loads.values()

    def get_transformers(self) -> ValuesView[Any]:
        """return a List-like view of network Transformer objects"""
        return self.transformers.values()

    def set_load_kw(self, load:str, kw:float) -> None:
        self.loads[load].set_kw(kw)

    def set_load_kvar(self, load:str, kvar:float) -> None:
        self.loads[load].set_kvar(kvar)

    def __str__(self) -> str:
        """return something informative"""
        df_dict = self.to_dataframes()
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
    def __init__(self, name: str = '') -> None:
        self.name = name
        self.phases = [False, False, False]
        self.parent = None # pointer to parent Node. only one parent for radial networks
        self.loads : List[Load] = []
        self.capacitors : List[Capacitor] = []
        self.sum_spu = np.zeros(3) # 3x1 complex array that holds the sum of all load.spu on this node
        self.sum_cappu = np.zeros(3) # 3x1 array that holds the sum of all capacitor.cappu on this node
        self.controller = None
        self.Sbase = 1000000.0
        self.Vbase = 0.0
        self.Ibase = 0.0
        self.Zbase = 0.0

    def __str__(self) -> str:
        return f"{self.name}, {self.phases}"

    def to_series(self) -> pd.Series:
        data = [self.name, self.phases, self.loads, self.controller]
        return pd.Series(data, self.series_index)


class Line:
    series_index = ['(tx,rx)', 'name', 'phases', 'config', 'length', 'FZpu']
    # TODO: might be helpful to include a list of pointers to all Lines in the class, and do the same for Node, etc.
    # see: http://effbot.org/pyfaq/how-do-i-get-a-list-of-all-instances-of-a-given-class.htm
    def __init__(self, network: Network, key: Tuple[str,str] = None, name: str = '') -> None:
        self.key = key # tuple of (txnode_name, rxnode_name)
        self.name = name # string, the name DSS uses to refer to this line
        self.phases = [False, False, False]
        self.config = None
        self.length = 0.0
        self.FZpu = np.zeros((3,3), dtype = 'complex')
        self.Sbase = 1000000.0
        self.Vbase = 0.0
        self.Ibase = 0.0
        self.Zbase = 0.0
        # add this Line to the network
        network.lines[key] = self
        tx, rx = key
        if rx not in network.adj[tx]:
            network.adj[tx].append(rx)  # add this as a Line in the adjacency list

    def __str__(self) -> str:
        return str(self.key)

    def to_series(self) -> pd.Series:
        data = [self.key, self.name, self.phases, self.config, self.length, self.FZpu]
        return pd.Series(data, self.series_index)


class Load:
    series_index = ['name', 'conn', 'phases', 'type', 'aPQ', 'aI', 'ppu', 'qpu', 'spu']
    def __init__(self, name: str = '') -> None:
        self.name = name
        self.conn = ''
        self.phases = [False, False, False]
        self.type = None
        self.aPQ = None
        self.aI = None
        self.aZ = None
        self.ppu = 0.0 + 0j
        self.qpu = 0.0
        self.spu = 0.0 + 0j
        self.kw = 0.0
        self.kvar = 0.0
        self.node : Union[ Node, Any ]= None

    def __str__(self) -> str:
        return self.name

    def set_kvar(self, kvar: float) -> None:
        """
        Reset a load's kvar and recalculate all load parameters based on new kvar.
        Note that this also updates the load's node's sum_spu
        """
        old_spu = self.spu
        self.kvar = kvar

        #divide by number of phases
        self.qpu = self.kvar / 1000 / self.phases.count(True)

        # set to 3x1 based on phases
        self.qpu = np.asarray([self.qpu if phase else 0 for phase in self.phases])
        self.spu = self.ppu + 1j*self.qpu

        # Update the sum_spu of the node that this load belongs to
        # Subtracting the old value from the sum, then add the current value
        self.node.sum_spu = np.subtract(self.node.sum_spu, old_spu)
        self.node.sum_spu = np.add(self.node.sum_spu, self.spu)

    def set_kw(self, kW: float) -> None:
        """
        Reset a load's kW and recalculate all load parameters based on new kvar.
        Note that this also updates the load's node's sum_spu
        """
        old_spu = self.spu
        self.kW = kW

        # divide by number of phases
        ppu = self.kW / 1000 / self.phases.count(True)

        # set to 3x1 based on phases
        self.ppu = np.asarray([ppu if phase else 0 for phase in self.phases])
        self.spu = self.ppu + 1j*self.qpu

        # Update the sum_spu of the node that this load belongs to
        # Subtracting the old value from the sum, then add the current value
        self.node.sum_spu = np.subtract(self.node.sum_spu, old_spu)
        self.node.sum_spu = np.add(self.node.sum_spu, self.spu)

    def to_series(self) -> pd.Series:
        data = [self.name, self.conn, self.phases, self.type, self.aPQ, self.aI, self.ppu, self.qpu, self.spu]
        return pd.Series(data, self.series_index)


class Controller:
    def __init__(self, name: str = '') -> None:
        self.node = None
        self.name = name
        self.phases = [False, False, False]
        self.wpu = np.zeroes(3)
        self.wmaxpu = np.zeroes(3)
        self.fes = np.zeroes(3)
        self.hpfes = np.zeroes(3)
        self.lpfes = np.zeroes(3)
        self.kintes = np.zeroes(3)


    def __str__(self) -> str:
        return self.name


class Capacitor:
    def __init__(self, name: str = '') -> None:
        self.name = name
        self.phases = [False, False, False]
        self.conn = ''
        self.cappu = np.zeros((3, 3), dtype='float')
    def __str__(self) -> str:
        return self.name


class Transformer (Line):
    def __init__(self, network: Network, key: Tuple[str, str], name: str, num_windings: int) -> None:
        super().__init__(network, key, name)  # initialize this as a Line with 0 length
        self.num_windings = num_windings
        self.Vbase = np.zeros((3,), dtype='complex')  # 3x1, initialize to nominal voltage
        self.kV = 0.0
        self.kVA = 0
        self.conn = ''  # 'wye' or 'delta

        # add this Transformer to the network
        network.transformers[name] = self

