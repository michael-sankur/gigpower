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
        self.iterations = -1  # stores number of iterations of FBS until convergence
        self.Vtest = np.zeros(3, dtype='complex')
        self.Vref = network.Vbase * np.ones(3, dtype='complex')
        self.solved_nodes = dict()
        self.solved_lines = dict()
        self.tolerance = -1  # stores the tolerance
        self.diff = -1  # stores the final value of Vtest - Vref at convergence

        """ Set up voltages and currents for all nodes """
        for node in network.get_nodes():
            node_dict = dict()
            self.solved_nodes[node.name] = node_dict
            # initialize a zeroed out array for each solution var
            node_dict['V'] = network.Vbase * np.ones(3, dtype='complex')
            node_dict['Inode'] = np.zeros(3, dtype='complex')
            node_dict['S'] = np.zeros(3, dtype='complex')
            node_dict['s'] = np.zeros(3, dtype = 'complex') # tracks voltage dependent load
        for line in network.get_lines():
            line_dict = dict()
            self.solved_lines[line.name] = line_dict
            line_dict['I'] = np.zeros(3, dtype='complex')

    def update_voltage_forward(network: Network, parent: Node, child: Node) -> None:
        """
        updates voltage at child based on parent according to:
        child_node_dict['V'] = parent_dict['V'] - (parent,child).FZpu * (parent,child).I
        """
        parent_dict = self.solved_nodes[parent.name]
        child_node_dict = self.solved_nodes[child_node.name]
        line_key = (parent.name, child_node.name)
        parent_V = parent_dict['V']
        child_V = child_node_dict['V']
        FZpu = network.lines[line_key].FZpu
        I = self.solved_lines[line_key]['I']

        # child node voltage = parent node voltage - current(child_node, parent)
        new_child_V = parent_V - np.matmul(FZpu, I)
        # zero out voltages for not existent phases at child node
        child_node_dict['V'] = mask_phases(new_child_V, child_node.phases)
        return None

    def update_voltage_backward(network: Network, child: Node, parent: Node) -> None:
        """
        updates voltage at parent node only for nodes existing on child node.
        """
        parent_dict = self.solved_nodes[parent.name]
        child_dict = self.solved_nodes[child.name]
        child_V = child_node_dict['V']
        line_key = (parent.name, child.name)
        FZpu = network.lines[line_key].FZpu
        I = self.solved_lines[line_key]['I']
        for phase_idx in range(3):
            if child.phases[phase_idx]: # if this phase is present on child
                parent_dict['V'][phase_idx] = child_V[phase_idx] + np.matmul(FZpu[phase_idx], I)
        # zero out voltages for not existent phases at parent node
        parent_dict['V'] = mask_phases(parent_dict['V'], parent.phases)
        return None

    def update_current(network: Network, line_in: Line) -> None:
        downstr_node_name = line_in.key[1]
        solved_line_dict = self.solved_lines[line_in.key]
        downstr_node_dict = self.solved_nodes[ downstr_node_name ]
        downstr_node_phases = network.nodes[ downstr_node_name ].phases
        line_I = solved_line_dict['I']
        line_s = solved_line_dict['s']
        downstr_node_V = downstr_node_dict['V']

        new_line_I = np.conj(np.divide(line_s, downstr_node_V))
        # TODO: confirm that np.divide is the same as matlab right divide './'
        new_line_I = mask_phases(new_line_I, downstr_node_phases)

    def update_voltage_dependent_load(self, network: Network) -> None:
        """
        update s at network loads
        """
        # s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V)).^2) - 1j * cappu + wpu
        aPQ = 1.00
        aI = 0
        aZ =  0
        # TODO; get aPQ, aI, aZ from dss file
        for node_name, node_dict in self.solved_nodes.items():
            node = network.node[node_name]
            node_V = node_dict['V']
            wpu = self.wpu = np.zeroes(3) # TODO: get wpu from dss file
            spu = node.load.spu
            node_dict['s'] = np.multiply(spu, aPQ + np.multiply(aI * abs(V)) ) + np.multiply(aZ (np.linalg.matrix_power(abs(V), 2))) - 1j * cappu + wpu

    def solved_nodes_df(self) -> Iterable:
        """
        returns solved nodes as a dataframe indexed by node
        """
        return pd.DataFrame.from_dict(self.solved_nodes).transpose()

    def solved_lines_df(self) -> Iterable:
        """
        returns solved lines as a dataframe indexed by line
        """
        return pd.DataFrame.from_dict(self.solved_lines).transpose()

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
        print(self.solved_nodes_df())
        print("Solved Lines:")
        print(self.solved_lines_df())
        print()
