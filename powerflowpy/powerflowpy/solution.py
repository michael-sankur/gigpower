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
        self.solved_nodes = dict()
        self.solved_lines = dict()
        self.tolerance = -1  # stores the tolerance
        self.diff = -1  # stores the final value of Vtest - Vref at convergence
        self.root = None # keeps a pointer to the root node

        """ Set up voltages and currents for all nodes """
        for node in network.get_nodes():
            node_dict = dict()
            self.solved_nodes[node.name] = node_dict
            # set initial voltage at each node to reference voltage
            node_dict['V'] = self.Vref
            # initialize a zeroed out array for each solution var
            node_dict['Inode'] = np.zeros(3, dtype='complex')
            node_dict['S'] = np.zeros(3, dtype='complex')
            node_dict['s'] = np.zeros(3, dtype = 'complex') # tracks voltage dependent load
        for line in network.get_lines():
            line_dict = dict()
            self.solved_lines[line.key] = line_dict
            line_dict['I'] = np.zeros(3, dtype='complex')

    def update_voltage_forward(self, network: Network, parent: Node, child: Node) -> None:
        """
        updates voltage at child based on parent according to:
        child_node_dict['V'] = parent_dict['V'] - (parent,child).FZpu * (parent,child).I
        """
        parent_dict = self.solved_nodes[parent.name]
        child_dict = self.solved_nodes[child.name]
        line_key = (parent.name, child.name)
        parent_V = parent_dict['V']
        child_V = child_dict['V']
        FZpu = network.lines[line_key].FZpu
        I = self.solved_lines[line_key]['I']

        # child node voltage = parent node voltage - current(child_node, parent)
        new_child_V = parent_V - np.matmul(FZpu, I)
        # zero out voltages for not existent phases at child node
        child_dict['V'] = mask_phases(new_child_V, (3,), child.phases)
        return None

    def update_voltage_backward(self, network: Network, child: Node) -> None:
        """
        updates voltage at parent node only for phases existing on child node.
        """
        parent = child.parent
        parent_dict = self.solved_nodes[parent.name]
        child_dict = self.solved_nodes[child.name]
        child_V = child_dict['V']
        line_key = (parent.name, child.name)
        FZpu = network.lines[line_key].FZpu
        I = self.solved_lines[line_key]['I']
        # update voltages at parent only for phases existing on child
        for phase_idx in range(3):
            if child.phases[phase_idx]: # if this phase is present on child
                parent_dict['V'][phase_idx] = child_V[phase_idx] + np.matmul(FZpu[phase_idx], I)
        # zero out voltages for not existent phases at parent node
        parent_dict['V'] = mask_phases(parent_dict['V'], (3,), parent.phases)
        return None

    def update_current(self, network: Network, line_in: Line) -> None:
        """
        Updates current on a line based on the downstream node. Used during backward sweep.
        """
        node_name = line_in.key[1]
        line_dict = self.solved_lines[line_in.key]
        node_dict = self.solved_nodes[ node_name ]
        node_phases = network.nodes[ node_name ].phases
        line_I = line_dict['I']
        node_s = node_dict['s']
        node_V = node_dict['V']
        new_line_I = np.conj(np.divide(node_s, node_V))
        # sum currents over all node's child segments
        for child_name in network.adj[node_name]:
            child_segment = (node_name, child_name)
            new_line_I = new_line_I + self.solved_lines[child_segment]['I']

        # TODO: confirm that np.divide is the same as matlab right divide './'
        new_line_I = mask_phases(new_line_I, (3,), node_phases)
        line_dict['I'] = new_line_I

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
        for node_name, node_dict in self.solved_nodes.items():
            node = network.nodes[node_name]
            node_V = node_dict['V']
            wpu = np.zeros(3) # TODO: get wpu from dss file
            cappu = np.zeros(3)  # TODO: get cappu from dss file
            spu = node.load.spu if node.load else np.zeros(3)
            node_dict['s'] = np.multiply(spu, aPQ + np.multiply(aI, abs(node_V)) ) + np.multiply(aZ, (np.power(abs(node_V), 2))) - 1j * cappu + wpu
        return None

    def solved_nodes_df(self) -> Iterable:
        """
        returns solved nodes as a dataframe indexed by node
        """
        return pd.DataFrame.from_dict(self.solved_nodes, orient = 'index')

    def solved_lines_df(self) -> Iterable:
        """
        returns solved lines as a dataframe indexed by line
        """
        return pd.DataFrame.from_dict(self.solved_lines, orient = 'index')

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
        print("Solved Nodes:")
        for key,d in self.solved_nodes.items():
            print(f"{key} \t {d['V']}")
        print("Solved Lines:")
        print(self.solved_lines_df())
        print()

