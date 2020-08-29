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

    def update_voltage(network: Network, source_node: Node, target_node: Node, direction: str) -> None:
        """
        updates voltage at target_node based on source_node according to:
        target_node_dict['V'] = source_node_dict['V'] - (source,target).FZpu * (source,target).I
        """
        source_node_dict = self.solved_nodes[source_node.name]
        target_node_dict = self.solved_nodes[target_node.name]
        line_key = (source_node.name, target_node.name)
        source_V = source_node_dict['V']
        target_V = target_node_dict['V']
        FZpu = network.lines[line_key].FZpu
        I = self.solved_lines[line_key]['I']
        # if forward, subtract Z*I, otherwise add Z*I
        sign = -1 if direction == 'forward' else 1

        # target node voltage = source node voltage +/- current(target_node, source_node)
        new_target_V = source_V + sign * np.matmul(FZpu, I)
        # zero out voltages for not existent phases at target node
        # ELAINE: move this out or toggle for forward/backward
        target_node_dict['V'] = mask_phases(new_target_V, target_node.phases)
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
