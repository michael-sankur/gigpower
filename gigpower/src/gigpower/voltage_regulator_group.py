from . line_group import LineGroup
from . circuit_element_group import CircuitElementGroup
from . voltage_regulator import VoltageRegulator, Tuple
import numpy as np
import copy
from typing import List

class VoltageRegulatorGroup(LineGroup):
    dss_module_name = 'RegControls'
    ele_class = VoltageRegulator

    def __init__(self, dss, bus_group):
        super().__init__(dss, bus_group=bus_group)
        related_vr = {}
        for n in range(len(dss.RegControls.AllNames())):
            dss.RegControls.Name(dss.RegControls.AllNames()[n])
            dss.Circuit.SetActiveBus(
                dss.CktElement.BusNames()[0].split(".")[0])
            if dss.Bus.Name() in related_vr.keys():
                related_vr[dss.Bus.Name()].append(n)
            else:
                related_vr[dss.Bus.Name()] = [n]
        self.voltage_regulator_index_dict = related_vr

    def add_element(self, vr):
        """
        Add the VR to the VRGroup, including all the lines of the VR 
        in the VRGroup's adjacency matrices
        """
        CircuitElementGroup.add_element(self, vr, inc_num_elements=True)
        try:
            self._key_to_element_dict[vr.key].append(vr)
        except KeyError:
            self._key_to_element_dict[vr.key] = [vr]
        super()._add_edge(vr.line_tx_to_reg)
        super()._add_edge(vr.line_reg_to_tx)

    def get_bus_ids(self, which: str):
        """
        overwrite super to get tx buses of upstream, downstream lines
        in vr index order
        param which = 'tx' or 'rx'
        resulting list length is 2 * self.num_elements
        even number indices are upstream_line tx bus indices
        odd number indices are downstream_line tx bus indices
        """
        buses = []
        for vr in self.get_elements():
            buses.append(getattr(vr.line_tx_to_reg, which))
            buses.append(getattr(vr.line_reg_to_tx, which))
        return buses

    def get_gain_matrix(self) -> np.ndarray:
        """
        n x 1 matrix of NEGATIVE gain values, indexed by voltage_regulator index
        Negative sign for voltage ratio purposes
        """
        return [-1 * vr.gamma for vr in self.get_elements()]

    def get_adj_set(self, reversed=False):
        """
        Return self.adj, with all voltage regulators between the same two buses
        into the same edge.
        ex: {'650':['rg60', 'rg60', 'rg60]} -> {'650': ['rg60']}
        """
        if reversed:
            new_adj = copy.deepcopy(self.reverse_adj)
        else:
            new_adj = copy.deepcopy(self.adj)
        for node, neighbors in new_adj.items():
            new_adj[node] = list(set(neighbors))
        return new_adj

    def get_idx(self, obj):
        """
        Handles tuple representation. If passed a tuple, returns a list of VRs
        on the Line represented by the tuple
        """
        if isinstance(obj, tuple):
            return self._key_to_element_dict[obj]  
        else:
            return CircuitElementGroup.get_idx(self, obj)
