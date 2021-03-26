from . line_group import LineGroup
from . circuit_element_group import CircuitElementGroup
from . voltage_regulator import VoltageRegulator
from typing import List
import numpy as np


class VoltageRegulatorGroup(LineGroup):
    dss_module_name = 'RegControls'
    ele_class = VoltageRegulator

    def __init__(self, dss, line_group):
        super().__init__(dss, bus_group=line_group.buses, line_group=line_group)
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
        super()._add_edge(vr.line_tx_to_reg)
        super()._add_edge(vr.line_reg_to_tx)

    # def get_bus_ids(self, which: str) -> List:
    #     """
    #     overwrite super to get tx buses of upstream, downstream lines
    #     in vr index order
    #     param which = 'tx' or 'rx'
    #     resulting list length is 2 * self.num_elements
    #     even number indices are upstream_line tx bus indices
    #     odd number indices are downstream_line tx bus indices
    #     """
    #     buses = []
    #     for vr in self.get_elements():
    #         buses.append(getattr(vr.upstream_line, which))
    #         buses.append(getattr(vr.downstream_line, which))
    #     return buses

    def get_gain_matrix(self) -> np.ndarray:
        """
        n x 1 matrix of NEGATIVE gain values, indexed by voltage_regulator index
        Negative sign for voltage ratio purposes
        """
        return [-1 * vr.gamma for vr in self.get_elements()]

        
