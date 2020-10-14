# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 28 August 2020
# Implement FBS to solve Power Flow for a radial distribution network.

from typing import Iterable, Dict, Tuple
from . utils import mask_phases
from . network import Network, Node, Line
import pandas as pd # type: ignore
import numpy as np # type: ignore

class Solution:

    def __init__(self, network: Network, tol: float = -1, max_iter: int = -1) -> None:

        # PARAMETERS
        self.iterations = 0  # stores number of iterations of FBS until convergence
        self.Vtest = np.zeros(3, dtype='complex')
        self.Vref = np.array( [1, np.exp(1j*240*np.pi/180), np.exp(1j*120*np.pi/180)], dtype = complex)
        self.tolerance = -1  # stores the tolerance at most recent completed iteration
        self.diff = -1  # stores the final value of Vtest - Vref at convergence
        # hard coded VREF = [1, 1*exp(j*240*pi/180), 1*exp(j*120*pi/180)]
        #TODO: set reference voltage based on source bus

        # SOLUTION VARIABLES
        self.V = dict() # 3x1 complex pu voltages phasors. 3 x num_nodes total datapoints. Indexed by node name.
        self.Inode = dict() # 3x1 complex pu current phasors delivered to node. 3 x num_nodes total data points. Indexed by node name.
        self.I = dict() # 3x1 complex pu current phasors. 3 x num_lines total datapoints. Indexed by line key.
        self.Stx : Dict[ Tuple[str,str], Iterable] = dict() # line transmitting end power. 3 x num_lines total datapoints. Indexed by line key.
        # 3x1 line transmitting end power. 3 x num_lines total datapoints. Indexed by line key.
        self.Srx : Dict[Tuple[str, str], Iterable] = dict()
        # Total power at each node. 3 x num_nodes total datapoints. Indexed by node name.
        self.sV: Dict[str, Iterable] = dict()

        # HELPER VARIABLES
        # self.s: intermediate values for voltage dependent loads. 3 x num_nodes total datapoints. Indexed by node name.
        self.s = dict()
        self.root = None # keeps a pointer to the root node
        self.network = network # keeps a pointer to the network that it is solving

        """ Set up voltages and currents for all nodes """
        for node in network.get_nodes():
            # set initial voltage at each node to reference voltage
            self.V[node.name] = np.zeros(3, dtype='complex')
            # initialize a zeroed out array for each solution var
            self.Inode[node.name] = np.zeros(3, dtype='complex')
            self.sV[node.name] = np.zeros(3, dtype='complex')
            # initialize voltage dependent load
            self.s[node.name] = np.zeros(3, dtype='complex')

        """ Set up currents for all lines """
        for line in network.get_lines():
            self.I[line.key] = np.zeros(3, dtype='complex')

        """ Initialize convergence tracking data frames """
        # TODO: setup index order based on topo order
        index = [ n.name for n in network.get_nodes() ]
        cols = [
            'curr_iter',
            'Vprev_a',
            'Vprev_b',
            'Vprev_c',
            'Vcurr_a',
            'Vcurr_b',
            'Vcurr_c',
            'max_absVdiff',
        ]
        self.V_forward_delta = pd.DataFrame( [], index, cols)
        self.V_backward_delta = pd.DataFrame([], index, cols)

    def update_voltage_forward(self, network: Network, parent: Node, child: Node) -> None:
        """
        updates voltage at child based on parent according to:
        child_node_dict['V'] = parent_dict['V'] - (parent,child).FZpu * (parent,child).I
        """
        line_key = (parent.name, child.name)

        # if line_key is a voltage_regulator:
        #     # voltage regulator a pair of lines: open with upstream bus, one with downstream bus
        #     # Vbase is the same on both sides
        #     # power going in = power going out
        #     child_V = gamma * parent_V # gamma is given
        #     child_V * downstream_line_I == parent_V * upstream_line_I
        #     for phase in phases:
        #         VR_downstream_line_voltage =  gamma * VR_upstream_voltage
        #         VR_downstream_line_current = VR_upstream_line_voltage * VR_upstream_line_current / VR_downstream_line_voltage
        #         # VR_downstream_line_voltage * VR_downstream_line_current = VR_upstream_line_voltage * VR_downstream_line_current

        parent_V = self.V[parent.name]
        FZpu = network.lines[line_key].FZpu
        I = self.I[line_key]

        # child node voltage = parent node voltage - current(child_node, parent)
        new_child_V = parent_V - np.matmul(FZpu, I)
        # zero out voltages for non-existant phases at child node
        self.V[ child.name ] = mask_phases(new_child_V, (3,), child.phases)
        #TODO: find out what -0j is. This happens after masking phases on the child.
        return None

    def update_voltage_backward(self, network: Network, child: Node) -> None:
        """
        updates voltage at parent node only for phases existing on child node.
        """
        parent = child.parent
        child_V = self.V[child.name]
        parent_V = self.V[parent.name]  # type: ignore
        line_key = (parent.name, child.name)  # type: ignore
        # if line_key is voltage_regulator:
        #     line2 = VR_downstream_line
        #     line1 = VR_upstream_line
        #     for phase in phases:
        #         VR_upstream_line_voltage = 1/ gamma * VR_downstream_voltage
        #         VR_upstream_line_current = VR_downstream_line_voltage * VR_downstream_line_current / VR_upstream_line_voltage
        #         # VR_downstream_line_voltage * VR_downstream_line_current = VR_upstream_line_voltage * VR_downstream_line_current


        FZpu = network.lines[line_key].FZpu
        I = self.I[line_key]
        # update voltages at parent only for phases existing on child
        for phase_idx in range(3):
            if child.phases[phase_idx]:  # if this phase is present on child
                parent_V[phase_idx] = child_V[phase_idx] + np.matmul(FZpu[phase_idx], I)
        # zero out voltages for non-existant phases at parent node
        self.V[parent.name] = mask_phases(parent_V, (3,), parent.phases) # type: ignore
        return None

    def update_parent_current(self, network: Network, line_in: Line) -> None:
        """
        Updates current at parent segment based on a line_in, according to:
        [parent_segment]---->parent_node---->[line_in]----->current_node
        Used during backward sweep
        """
        parent_name = line_in.key[0] # type: ignore
        parent_node = self.network.nodes.get(parent_name)
        parent_seg_key = (parent_node.parent.name, parent_name) # type: ignore
        parent_seg_phases = network.lines[parent_seg_key].phases
        line_in_I = self.I[line_in.key]
        parent_s = self.s[parent_name]
        parent_V = self.V[parent_name]
        # np.divide produces a NaN for positions at which node_V is 0 because the phases are not existant on node
        new_parent_I = line_in_I + np.conj(np.divide(parent_s, parent_V))
        new_parent_I = mask_phases(new_parent_I, (3,), parent_seg_phases)
        self.I[parent_seg_key] = new_parent_I

    def update_current(self, network: Network, line_in: Line) -> None:
        """
        Updates current on a line_in based on the downstream node, according to:
        upstream_node--->[line_in]---> downstream_node ===> [0 or many lines out] ===> 0 or many child nodes
        Used during backward sweep.
        """
        node_name = line_in.key[1] # type: ignore
        line_phases = network.lines[ line_in.key ].phases #type: ignore
        node_s = self.s[node_name]
        node_V = self.V[node_name]
        # np.divide produces a NaN for positions at which node_V is 0 because the phases are not existant on node
        new_line_I = np.conj(np.divide(node_s, node_V)) # TODO: consider dividing manually only by phases on the node
        # sum currents over all node's child segments
        for child_name in network.adj[node_name]:
            child_segment = (node_name, child_name)
            new_line_I = new_line_I + self.I[child_segment]

        new_line_I = mask_phases(new_line_I, (3,), line_phases)
        self.I[line_in.key] = new_line_I

    def update_voltage_dependent_load(self, node: Node) -> None:
        """
        updates voltage_dependent_load (Solution.s[node]) at the input node.
        Used during backward sweep.
        """
        # s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V)).^2) - 1j * cappu + wpu
        aPQ = np.ones(3)
        aI = np.zeros(3)
        aZ = np.zeros(3)
        # TODO; get aPQ, aI, aZ from dss file
        # dss.LoadModels - above equal to constant p and q. We use model 1 and 8

        node_V = self.V[node.name]
        wpu = np.zeros(3) # TODO: will be set as argument
        cappu = node.sum_cappu
        spu = node.sum_spu
        self.s[node.name] = np.multiply(spu, aPQ + np.multiply(aI, abs(node_V)) ) + np.multiply(aZ, (np.power(abs(node_V), 2))) - 1j * cappu + wpu
        return None

    def calc_S(self) -> None:
        """ Calculate Stx and Srx """
        for line in self.network.get_lines():
            tx_name, rx_name = line.key
            self.Stx [ line.key ] = np.multiply(self.V[ tx_name ], np.conj(self.I[line.key]))
            self.Srx [ line.key ] = np.multiply(self.V[ rx_name ], np.conj(self.I[line.key]))

    def calc_sV(self) -> None:
        """ Final calculation of voltage dependent complex loads. """
        # TODO: # handle multiple loads with update_s method. This is redundant (equivalent to update s)
        for node in self.network.get_nodes():
            self.update_voltage_dependent_load(node) # update self.s one last time
        self.sV = self.s  # set self.sv to self.s

    def calc_Inode(self) -> None:
        """ Calculate self.Inode (currents consumed at each node) """
        for node in self.network.get_nodes():
            node_V = self.V[ node.name ]
            node_sV = self.sV[ node.name ]
            node_I = np.conj(np.divide(node_sV,node_V))
            self.Inode[ node.name ] = mask_phases(node_I, (3,), node.phases)

    def V_df(self) -> Iterable:
        """
        returns self.V as a dataframe indexed by node name
        """
        return pd.DataFrame.from_dict(self.V, orient = 'index', columns=['A', 'B', 'C'])

    def I_df(self) -> Iterable:
        """
        returns self.I as a dataframe indexed by line name
        """
        Idf = pd.DataFrame.from_dict(self.I, orient='index', columns=['A', 'B', 'C'])
        # reindex lines to match opendss file
        new_index = ( [ self.network.lines.get(k).name for k in self.I.keys()] ) # type: ignore
        Idf.index = new_index
        return Idf

    def Inode_df(self) -> Iterable:
        """
        returns self.Inode as a dataframe indexed by node name
        """
        return pd.DataFrame.from_dict(self.Inode, orient='index', columns=['A', 'B', 'C'])

    def Stx_df(self) -> Iterable:
        """
        returns self.Stx as a dataframe indexed by line name
        """
        Stx = pd.DataFrame.from_dict(
            self.Stx, orient='index', columns=['A', 'B', 'C'])
        # reindex lines to match opendss file
        new_index = ([self.network.lines.get(k).name for k in self.I.keys()]) # type: ignore
        Stx.index = new_index
        return Stx

    def Srx_df(self) -> Iterable:
        """
        returns self.Srx as a dataframe indexed by line name
        """
        Srx = pd.DataFrame.from_dict(
            self.Srx, orient='index', columns=['A', 'B', 'C'])
        # reindex lines to match opendss file
        new_index = ([self.network.lines.get(k).name for k in self.I.keys()]) # type: ignore
        Srx.index = new_index
        return Srx

    def sV_df(self) -> Iterable:
        """
        returns self.sV as a dataframe indexed by node name
        """
        return pd.DataFrame.from_dict(self.sV, orient='index', columns=['A', 'B', 'C'])

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
        print("\n Parameters:")
        print(self.params_df())

        print("\n V solution")
        print(self.V_df())

        print("\n I solution")
        print(self.I_df())

        print("\n Inode solution")
        print(self.Inode_df())

        print("\n Stx solution")
        print(self.Stx_df())

        print("\n Srx solution")
        print(self.Srx_df())

        print("\n sV solution")
        print(self.sV_df())
        print()

