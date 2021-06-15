# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: March 21 2021
# Create FBS Solution class, a namespace for calculations used by fbs

import numpy as np  # type: ignore
from . voltage_regulator import VoltageRegulator
from typing import List, Dict, Iterable
from . solution import Solution
from . utils import mask_phases, topo_sort, get_reverse_adj
import copy


class SolutionFBS(Solution):
    
    @classmethod
    def set_zip_values(cls, zip_v):
        """
        sets zip values for the Solution class
        param zip_V: List or nd.array with 7 values
        [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min voltage pu]
        Note that zip values are set both on the Solution class and Circuit
        class
        """
        Solution.set_zip_values(zip_v)    

    def __init__(self, dss_fp: str, **kwargs):
        super().__init__(dss_fp, **kwargs)  # sets self.circuit
        self.set_adj()
        self.topo_order = topo_sort(self.circuit.buses.all_names(),
                                    self.adj)
        # save a pointer to the root bus
        self.root = self.circuit.buses.get_element(self.topo_order[0])
        self.Vref = np.array(
            [1, np.exp(1j*240*np.pi/180), np.exp(1j*120*np.pi/180)], dtype=complex)
        # TODO: move into class
        self.tolerance = abs((self.Vref[1]) * 10**-6)

    def set_adj(self):
        """
        Sets self.adj and self.reverse_adj for FBS algorithm
        FBS considers transformers and voltage regulators as edges, ignoring
        voltage regulator connections going upstream
        (otherwise, the circuit has cycles and is non-radial)
        """
        # voltage regulators between the same pair of buses
        # are consolidated into one edge by VoltageRegulatorGroup.get_adj_set()
        self.adj = {**self.circuit.lines.adj} # shallow copy of lines adjacency

        for node, xfm_adj in self.circuit.transformers.adj.items():
            # deep copy of adjacency lists for that need to include transformers
            if node not in self.adj: 
                self.adj[node] = xfm_adj
            else:
                self.adj[node] = copy.deepcopy(self.adj[node])
                self.adj[node] += xfm_adj
        vr_set_edges = self.circuit.voltage_regulators.get_adj_set()

        for node, vr_adj in vr_set_edges.items():
            # deep copy of adjacency lists for that need to include vrs
            if node not in self.adj:
                self.adj[node] = vr_adj
            else:
                self.adj[node] = copy.deepcopy(self.adj[node])
                self.adj[node] += vr_adj

        # remove any VR edges going upstream
        for vr in self.circuit.voltage_regulators.get_elements():
            tx, reg = vr.key
            try:
                self.adj[reg].remove(tx)
            except Exception:
                pass
        self.reverse_adj = get_reverse_adj(self.adj)

    def solve(self):
        """
        solve powerflow for self.circuit,
        save solution as self.V = bus complex voltage, self.I = line currents,
        self.S = bus powers
        based on fbs/fbs/fbs.py written by @elaguerta
        """
        topo_order = self.topo_order
        root_idx = self.circuit.buses.get_idx(self.root)
        converged = max(abs(self.Vtest - self.Vref)) <= self.tolerance

        while not converged and self.iterations < self.maxiter:
            # set V.root to Vref
            self.V[root_idx] = np.copy(self.Vref)  # type: ignore

            # FORWARD SWEEP: for bus in topo_order:
            for bus_name in topo_order:
                bus = self.circuit.buses.get_element(bus_name)
                if bus_name in self.adj:
                    children = self.adj[bus_name]
                    for child_name in children:
                        child = self.circuit.buses.get_element(child_name)
                        self.update_voltage_forward(bus, child)

            # update s at all buses
            self.calc_sV()

            # BACKWARD SWEEP: for node in reverse topo_order:
            for bus_name in reversed(topo_order):
                bus = self.circuit.buses.get_element(bus_name)
                # if this is a terminal node or junction node (not the root)
                parents = self.reverse_adj[bus_name]
                if len(parents) > 1:
                    raise ValueError(f"Network not radial. Bus {bus_name} has \
                        parents {parents}")
                elif len(parents) == 1:
                    parent = self.circuit.buses.get_element(parents[0])
                    # update sV at this bus, and at parent
                    self.calc_sV(bus)
                    self.calc_sV(parent)
                    line_in = self.circuit.lines.get_element(
                        (parent.__name__, bus.__name__))
                    # update line_in segment
                    self.update_current(line_in)
                    # update voltage at parent
                    self.update_voltage_backward(parent, bus)

            self.iterations += 1
            # set Vtest to the root's voltage
            self.Vtest = self.V[root_idx]  # type: ignore
            self.convergence_diff = max(abs(self.Vtest - self.Vref))
            # check convergence
            converged = self.convergence_diff <= self.tolerance

        # after convergence, update sV with V values at convergence 
        # and final calculations
        self.calc_sV()
        self.calc_Vmag()
        self.calc_Srx()
        self.calc_Stx()
        self.converged = converged

    def vr_forward(self, vr_list):
        """
        Enforce voltage regulator equations by phase.
        Voltage ratio equation: tx_node.V = gamma @ reg_node.V
        Conservation of power:
        reg_node.V @ [current entering reg_node]* ==
        tx_node.V @ [current entering tx_node]*
        [current @ reg_node]* =
        tx_node.V @ [current entering tx_node]*reg_node.V
        """
        # get all voltage regulators on the same edge
        for vr in vr_list:
            tx, reg = vr.key
            tx_idx = self.circuit.buses.get_idx(tx)
            reg_idx = self.circuit.buses.get_idx(reg)
            phases = vr.get_ph_idx_matrix()
            gamma = vr.gamma
            tx_V = self.V[tx_idx]
            reg_V = self.V[reg_idx]
            Itx = vr.Itx  # current at/leaving the tx node
            Ireg = vr.Ireg  # current entering reg_node
            reg_V[phases] = gamma * tx_V[phases]
            Ireg[phases] = 1/gamma * Itx[phases]
   
    def update_voltage_forward(self, parent, child):
        """
        updates voltage at child based on parent according to:
        child_node_dict['V'] = parent_dict['V'] - (parent,child).FZpu * (parent,child).I
        """
        line_key = (parent.__name__, child.__name__)
        line_out = self.circuit.lines.get_element(line_key)

        if isinstance(line_out, list):
            self.vr_forward(line_out)

        else:
            parent_idx = self.circuit.buses.get_idx(parent)
            child_idx = self.circuit.buses.get_idx(child)
            line_out_idx = self.circuit.lines.get_idx(line_key)
            parent_V = self.V[parent_idx]
            # child voltage = parent voltage - current(child, parent)
            FZpu = line_out.FZpu
            I_forward = self.I[line_out_idx]
            new_child_V = parent_V - np.matmul(FZpu, I_forward)
            # update V at child
            # zero out voltages for non-existant phases at child node
            self.V[child_idx] = new_child_V * child.phase_matrix

    def vr_backward(self, vr_list):
        """
        Enforce voltage regulator equations by phase.
        Voltage ratio equation: reg_node.V = 1/gamma @ tx_node.V
        Conservation of power:
        reg_node.V @ [current entering reg_node]* ==
        tx_node.V @ [current entering tx_node]*
        [current @ tx_node]* =
        reg_node.V @ [current entering reg_node]* / tx_node.V
        """
        # get all voltage regulators on the same edge
        for vr in vr_list:
            tx, reg = vr.key
            tx_idx = self.circuit.buses.get_idx(tx)
            reg_idx = self.circuit.buses.get_idx(reg)
            phases = vr.get_ph_idx_matrix()
            gamma = vr.gamma
            tx_V = self.V[tx_idx]
            reg_V = self.V[reg_idx]
            tx_V[phases] = 1/gamma * reg_V[phases]

    def update_voltage_backward(self, parent, child):
        """
        updates voltage at parent only for phases existing on child.
        """
        child_idx = self.circuit.buses.get_idx(child)
        parent_idx = self.circuit.buses.get_idx(parent)
        child_V = self.V[child_idx]
        line_key = (parent.__name__, child.__name__)  # type: ignore
        line_in = self.circuit.lines.get_element(line_key)
        
        if isinstance(line_in, list):
            self.vr_backward(line_in)

        else:
            line_in_idx = self.circuit.lines.get_idx(line_in)
            parent_V = self.V[parent_idx]
            FZpu = line_in.FZpu
            I_backwards = self.I[line_in_idx]
            # update voltages at parent only for phases existing on child
            phases = child.get_ph_idx_matrix()
            parent_V[phases] = child_V[phases] + np.matmul(FZpu[phases], I_backwards)
            # for phase_idx in range(3):
            #     if child.phase_matrix[phase_idx] == 1:  # if this phase is present on child
            #         parent_V[phase_idx] = child_V[phase_idx] + np.matmul(FZpu[phase_idx], I_backwards)
            # zero out non-existant phases at parent, save parent V
            self.V[parent_idx] = parent_V * parent.phase_matrix

    def vr_current(self, vr_list):
        """
        Enforce voltage regulator equations by phase.
        Voltage ratio equation: reg_node.V = 1/gamma @ tx_node.V
        Conservation of power: reg_node.V @ [current entering reg_node]* ==
        tx_node.V @ [current entering tx_node]*
        [current @ tx_node]* =
        reg_node.V @ [current entering reg_node]* / tx_node.V
        """
        # update all vrs on the same edge
        for vr in vr_list:
            line_in_idx = self.circuit.lines.get_idx(vr)
            tx, reg = vr.key
            phases = vr.get_ph_idx_matrix()
            vr.Ireg = self.I[line_in_idx]

            Itx = vr.Itx  # current entering/leaving the tx node
            Ireg = vr.Ireg  # current entering/leaving the reg_node
            # gamma = 1
            gamma = vr.gamma
            Itx[phases] = gamma * Ireg[phases]

    def update_current(self, line_in):
        """
        Updates current on a line_in based on the child node, according to
        tx_bus -[line_in]-> rx_bus =[0 or many lines out]=> 0 or many child buses
        Used during backward sweep.
        """
        if isinstance(line_in, list):
            self.vr_current(line_in)
        else:
            rx_bus_name = line_in.key[1]  # type: ignore
            rx_bus_idx = self.circuit.buses.get_idx(rx_bus_name)
            line_in_idx = self.circuit.lines.get_idx(line_in)
            # divide sV by V only where V is non-zero
            new_line_I = np.zeros(3, dtype=complex)
            non_zero_idx = np.where(self.V[rx_bus_idx] != 0)
            new_line_I[non_zero_idx] = np.conj(
                np.divide(self.sV[rx_bus_idx, non_zero_idx],
                          self.V[rx_bus_idx, non_zero_idx]))
            # sum currents over all node's child segments
            if rx_bus_name in self.adj:
                for child_name in self.adj[rx_bus_name]:
                    line_out = self.circuit.lines.get_element((rx_bus_name, child_name))
                    if isinstance(line_out, list):
                        # update I using all vrs on the same outgoing edge
                        for vr in line_out:
                            new_line_I += vr.Itx
                    else:  # update I summing over all outgoing lines
                        line_out_idx = self.circuit.lines.get_idx(line_out)
                        new_line_I = new_line_I + self.I[line_out_idx]
            # zero out non-existant phases
            # TODO: should SyntheticLines have phases? the old FBS would zero out I for
            # all SyntheticLines
            self.I[line_in_idx] = np.nan_to_num(new_line_I * line_in.phase_matrix)
