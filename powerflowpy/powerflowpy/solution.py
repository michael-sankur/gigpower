# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 28 August 2020
# Implement FBS to solve Power Flow for a radial distribution network.

from typing import Iterable, List, Dict
from . utils import mask_phases
from . network import *
import pandas as pd

class Solution:

    def __init__(self, network: Network, tol: float = -1, max_iter: int = -1) -> None:
        self.iterations = 0  # stores number of iterations of FBS until convergence
        self.Vtest = np.zeros(3, dtype='complex')
        self.Vref = np.array( [1, np.exp(1j*240*np.pi/180), np.exp(1j*120*np.pi/180)], dtype = complex)
        # hard coded VREF = [1, 1*exp(j*240*pi/180), 1*exp(j*120*pi/180)]
        #TODO: set reference voltage based on source bus
        self.V = dict() # 3x1 complex pu voltages phasors. Indexed by node name.
        self.I = dict() # 3x1 complex pu current phasors. Indexed by node name.
        self.Inode = dict() # 3x1 complex pu current phasors delivered to node. Indexed by node name.
        self.I = dict() # 3x1 complex pu current phasors. Indexed by line key.
        self.S = dict() # 3x1 complex pu power phasors delivered to node. Indexed by node name.
        self.s = dict() # 3x1 stores intermediate values for voltage dependent loads
        self.sV = dict() # 3x1 complex pu power phasors consumed by node. Indexed by node name.
        self.tolerance = -1  # stores the tolerance
        self.diff = -1  # stores the final value of Vtest - Vref at convergence
        self.root = None # keeps a pointer to the root node

        """ Set up voltages and currents for all nodes """
        for node in network.get_nodes():
            # set initial voltage at each node to reference voltage
            self.V[node.name] = self.Vref
            # initialize a zeroed out array for each solution var
            self.Inode[node.name] = np.zeros(3, dtype='complex')
            self.S[node.name] = np.zeros(3, dtype='complex')
            self.sV[node.name] = np.zeros(3, dtype='complex')
            self.s[node.name] = np.zeros(3, dtype='complex')

        """ Set up currents for all lines """
        for line in network.get_lines():
            self.I[line.key] = np.zeros(3, dtype='complex')

    def update_voltage_forward(self, network: Network, parent: Node, child: Node) -> None:
        """
        updates voltage at child based on parent according to:
        child_node_dict['V'] = parent_dict['V'] - (parent,child).FZpu * (parent,child).I
        """
        line_key = (parent.name, child.name)
        parent_V = self.V[parent.name]
        child_V = self.V[child.name]
        FZpu = network.lines[line_key].FZpu
        I = self.I[line_key]

        # child node voltage = parent node voltage - current(child_node, parent)
        new_child_V = parent_V - np.matmul(FZpu, I)
        # zero out voltages for not existent phases at child node
        self.V[ child.name ] = mask_phases(new_child_V, (3,), child.phases)
        return None

    def update_voltage_backward(self, network: Network, child: Node) -> None:
        """
        updates voltage at parent node only for phases existing on child node.
        """
        parent = child.parent
        child_V = self.V[child.name]
        parent_V = self.V[parent.name]
        line_key = (parent.name, child.name)
        FZpu = network.lines[line_key].FZpu
        I = self.I[line_key]
        # update voltages at parent only for phases existing on child
        for phase_idx in range(3):
            if child.phases[phase_idx]: # if this phase is present on child
                parent_V[phase_idx] = child_V[phase_idx] + np.matmul(FZpu[phase_idx], I)
        # zero out voltages for not existent phases at parent node
        self.V[parent.name] = mask_phases(parent_V, (3,), parent.phases)
        return None

    def update_current(self, network: Network, line_in: Line) -> None:
        """
        Updates current on a line based on the downstream node. Used during backward sweep.
        """
        node_name = line_in.key[1]
        node_phases = network.nodes[ node_name ].phases
        line_I = self.I[line_in.key]
        node_s = self.s[node_name]
        node_V = self.V[node_name]
        new_line_I = np.conj(np.divide(node_s, node_V))
        # sum currents over all node's child segments
        for child_name in network.adj[node_name]:
            child_segment = (node_name, child_name)
            new_line_I = new_line_I + self.I[child_segment]

        # TODO: confirm that np.divide is the same as matlab right divide './'
        new_line_I = mask_phases(new_line_I, (3,), node_phases)
        self.I[line_in.key] = new_line_I

    def update_voltage_dependent_load(self, network: Network) -> None:
        """
        update s at network loads
        """
        #TODO: Recalculate s all at once over a supermatrix?

        # s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V)).^2) - 1j * cappu + wpu
        aPQ = 1.00
        aI = 0
        aZ =  0
        # TODO; get aPQ, aI, aZ from dss file
        for node in network.get_nodes():
            node_V = self.V[node.name]
            wpu = np.zeros(3) # TODO: get wpu from dss file
            cappu = np.zeros(3)  # TODO: get cappu from dss file
            spu = node.load.spu if node.load else np.zeros(3)
            self.S[node.name] = np.multiply(spu, aPQ + np.multiply(aI, abs(node_V)) ) + np.multiply(aZ, (np.power(abs(node_V), 2))) - 1j * cappu + wpu
        return None

    def V_df(self) -> Iterable:
        """
        returns self.V as a dataframe indexed by node name
        """
        return pd.DataFrame.from_dict(self.V, orient = 'index')

    def I_df(self) -> Iterable:
        """
        returns self.I as a dataframe indexed by line key
        """
        return pd.DataFrame.from_dict(self.I, orient = 'index')

    def Inode_df(self) -> Iterable:
        """
        returns self.Inode as a dataframe indexed by node name
        """
        return pd.DataFrame.from_dict(self.Inode, orient='index')

    def S_df(self) -> Iterable:
        """
        returns self.S as a dataframe indexed by node name
        """
        return pd.DataFrame.from_dict(self.S, orient='index')

    def sV_df(self) -> Iterable:
        """
        returns self.sV as a dataframe indexed by node name
        """
        return pd.DataFrame.from_dict(self.sV, orient='index')

    def params_df(self) -> Iterable:
        """
        returns solution paramaters as a dataframe
        """
        index = ['iterations', 'Vtest', 'Vref', 'tolerance', 'diff']
        data = [self.iterations, self.Vtest, self.Vref, self.tolerance, self.diff]
        return pd.DataFrame(data, index).transpose()

    def print_solution(self) -> None:
        """
        prints solution to stdout
        """
        print("Parameters:")
        print(self.params_df())

        print("V solution")
        print(self.V_df())

        print("I solution")
        print(self.I_df())

        print("Inode solution")
        print(self.Inode_df())

        print("S solution")
        print(self.S_df())

        print("sV solution")
        print(self.sV_df())
        print()

