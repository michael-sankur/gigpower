# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: March 21 2021
# Create FBS Solution class, a namespace for calculations used by fbs

import numpy as np  # type: ignore
from typing import List, Dict, Iterable
from . solution import Solution
from . utils import mask_phases, calculate_sV, topo_sort, get_reverse_adj
import copy

class SolutionFBS(Solution):

    def __init__(self, dss_fp: str):
        super().__init__(dss_fp)  # sets self.circuit
        # clean up voltage regulator topology
        self.adj = self._clean_adj(self.circuit.lines.adj, self.circuit.voltage_regulators)
        self.reverse_adj = get_reverse_adj(self.adj)

        self.topo_order = topo_sort(self.circuit.buses.all_names(),
                                    self.adj)
        # save a pointer to the root bus
        self.root = self.circuit.buses.get_element(self.topo_order[0])
        self.Vref = np.array(
            [1, np.exp(1j*240*np.pi/180), np.exp(1j*120*np.pi/180)], dtype=complex)
        # TODO: move into class
        self.tolerance = abs((self.Vref[1]) * 10**-6)

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

        # ask self.circuit to calculate these once
        self.spu = self.circuit.get_spu_matrix().transpose()
        self.aPQ = self.circuit.get_aPQ_matrix().transpose()
        self.aI = self.circuit.get_aI_matrix().transpose()
        self.aZ = self.circuit.get_aZ_matrix().transpose()
        self.wpu = self.circuit.get_wpu_matrix().transpose()
        self.cappu = self.circuit.get_cappu_matrix().transpose()
 
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
            self.update_sV()

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
                    self.update_sV(bus)
                    self.update_sV(parent)
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
        self.update_sV()
        self.converged = converged

    def update_voltage_forward(self, parent, child):
        """
        updates voltage at child based on parent according to:
        child_node_dict['V'] = parent_dict['V'] - (parent,child).FZpu * (parent,child).I
        """
        line_key = (parent.__name__, child.__name__)
        line_out = self.circuit.lines.get_element(line_key)

        if isinstance(line_out, List):  # list of Synthetic Lines indicates VRs
            for vr_line in line_out:
                for vr in vr_line.voltage_regulators:
                    tx, reg = vr.key
                    tx_idx = self.circuit.buses.get_idx(tx)
                    reg_idx = self.circuit.buses.get_idx(reg)
                    phases = vr.get_ph_idx_matrix()
                    gamma = vr.gamma
                    tx_V = self.V[tx_idx]
                    reg_V = self.V[reg_idx]
                    Itx = vr.Itx  # current at/leaving the tx node
                    Ireg = vr.Ireg  # current entering reg_node
                    # Enforce voltage regulator equations by phase.
                    # Voltage ratio equation: tx_node.V = gamma @ reg_node.V
                    # Conservation of power: reg_node.V @ [current entering reg_node]* == tx_node.V @ [current entering tx_node]*
                    # [current @ reg_node]* = tx_node.V @ [current entering tx_node]* / reg_node.V
                    reg_V[phases] = gamma * tx_V[phases]
                    Ireg[phases] = 1/gamma * Itx[phases]

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

    def update_voltage_backward(self, parent, child):
        """
        updates voltage at parent only for phases existing on child.
        """
        child_idx = self.circuit.buses.get_idx(child)
        parent_idx = self.circuit.buses.get_idx(parent)
        child_V = self.V[child_idx]
        line_key = (parent.__name__, child.__name__)  # type: ignore
        line_in = self.circuit.lines.get_element(line_key)
        
        if isinstance(line_in, list):  # implies voltage regulators
            for vr_line in line_in:
                line_in_idx = self.circuit.lines.get_idx(vr_line)
                for vr in vr_line.voltage_regulators:
                    tx, reg = vr.key
                    tx_idx = self.circuit.buses.get_idx(tx)
                    reg_idx = self.circuit.buses.get_idx(reg)
                    phases = vr.get_ph_idx_matrix()
                    gamma = vr.gamma
                    tx_V = self.V[tx_idx]
                    reg_V = self.V[reg_idx]
                    # Enforce voltage regulator equations by phase.
                    # Voltage ratio equation: reg_node.V = 1/gamma @ tx_node.V
                    # Conservation of power: reg_node.V @ [current entering reg_node]* == tx_node.V @ [current entering tx_node]*
                    # [current @ tx_node]* = reg_node.V @ [current entering reg_node]* / tx_node.V
                    #  update tx voltage per phase
                    tx_V[phases] = 1/gamma * reg_V[phases]

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

    def update_current(self, line_in):
        """
        Updates current on a line_in based on the child node, according to
        tx_bus -[line_in]-> rx_bus =[0 or many lines out]=> 0 or many chidl buses
        Used during backward sweep.
        """

        if isinstance(line_in, list):
            for vr_line in line_in:
                line_in_idx = self.circuit.lines.get_idx(vr_line)
                for vr in vr_line.voltage_regulators:
                    tx, reg = vr.key
                    phases = vr.get_ph_idx_matrix()
                    vr.Ireg = self.I[line_in_idx]

                    Itx = vr.Itx  # current entering/leaving the tx node
                    Ireg = vr.Ireg  # current entering/leaving the reg_node
                    # gamma = 1
                    gamma = vr.gamma
                    # Enforce voltage regulator equations by phase.
                    # Voltage ratio equation: reg_node.V = 1/gamma @ tx_node.V
                    # Conservation of power: reg_node.V @ [current entering reg_node]* == tx_node.V @ [current entering tx_node]*
                    # [current @ tx_node]* = reg_node.V @ [current entering reg_node]* / tx_node.V
                    Itx[phases] = gamma * Ireg[phases]
        else:
            rx_bus_name = line_in.key[1]  # type: ignore
            rx_bus_idx = self.circuit.buses.get_idx(rx_bus_name)
            line_in_idx = self.circuit.lines.get_idx(line_in)
            new_line_I = np.conj(np.divide(self.sV[rx_bus_idx], self.V[rx_bus_idx]))
            
            # sum currents over all node's child segments
            if rx_bus_name in self.adj:
                for child_name in self.adj[rx_bus_name]:
                    line_out = self.circuit.lines.get_element((rx_bus_name, child_name))
                    if isinstance(line_out, list):
                        for vr_line in line_out:
                            line_out_idx = self.circuit.lines.get_idx(vr_line)
                            for vr in vr_line.voltage_regulators:
                                new_line_I += vr.Itx
                    else:
                        line_out_idx = self.circuit.lines.get_idx(line_out)
                        new_line_I = new_line_I + self.I[line_out_idx]
            # zero out non-existant phases
            # TODO: should SyntheticLines have phases? the old FBS would zero out I for
            # all SyntheticLines
            self.I[line_in_idx] = np.nan_to_num(new_line_I * line_in.phase_matrix)

    def update_sV(self, bus=None):
        """
        TODO: this may need to move up to Solution
        updates self.sV (voltage dependent loads, or powers at each bus) 
        based on self.V's present values
        param bus: if Bus is given, will update self.sV only at the bus
        otherwise, updates all buses
        """
        # s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V)).^2) - 1j * cappu + wpu
        if not bus:  # all buses
            spu = self.spu
            aPQ, aI, aZ = self.aPQ, self.aI, self.aZ
            V, cappu, wpu = self.V, self.cappu, self.wpu
            self.sV = calculate_sV(V, spu, aPQ, aI, aZ, cappu, wpu)
        else:  # update at the given bus only
            bus_idx = self.circuit.buses.get_idx(bus)
            spu = self.spu[bus_idx]
            aPQ, aI, aZ = self.aPQ[bus_idx], self.aI[bus_idx], self.aZ[bus_idx]
            V, cappu, wpu = self.V[bus_idx], self.cappu[bus_idx], self.wpu[bus_idx]
            self.sV[bus_idx] = calculate_sV(V, spu, aPQ, aI, aZ, cappu, wpu)

    def _clean_adj(self, old_adj, vr_group):
        """
        cleans voltage regulators in an adjacency list,
        by representing all voltage regulators on the same line as one edge
        going downstream only
        (otherwise, the circuit has cycles and is non-radial)
        """
        new_adj = copy.deepcopy(old_adj)
        # remove all voltage regulator edges
        for vr in vr_group.get_elements():
            tx, reg = vr.key
            new_adj[tx].remove(reg)
            new_adj[reg].remove(tx)
        # add back edges going upstream only, condensed
        for vr in vr_group.get_elements():
            tx, reg = vr.key
            if reg not in new_adj[tx]:
                new_adj[tx].append(reg)

        return new_adj
