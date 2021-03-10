# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Circuit class to mirror a dss Circuit object
# used by Solution objects to solve powerflow

import numpy as np

from bus_group import BusGroup
from capacitor_group import CapacitorGroup
from line_group import LineGroup
from load_group import LoadGroup
from transformer_group import TransformerGroup
from voltage_regulator_group import VoltageRegulatorGroup
from utils import parse_phase_matrix


class Circuit():

    def __init__(self, dss, Sbase=10**6):
        """ initialize Circuit from an opendss object's current state"""
        self.Sbase = Sbase
        self.buses = BusGroup(dss)
        self.lines = LineGroup(dss)
        self.loads = LoadGroup(dss)
        self.capacitors = CapacitorGroup(dss)
        self.voltage_regulators = VoltageRegulatorGroup(dss, self.lines)
        self.transformers = TransformerGroup(dss, self.lines)

        self._assign_to_buses(self.loads)
        self._assign_to_buses(self.capacitors)
        self._assign_to_buses(self.voltage_regulators)
        self._assign_to_buses(self.transformers)

    def set_kW(self, load_name: str, kW: float):
        """
        sets a new kW for the given Load.
        Updates Load.spu, Load.ppu, Load.qpu, and Bus.sum_spu
        """
        load = self.loads.get_ckt_element(load_name)
        bus = self.bus.get_ckt_element(load.bus_name)
        old_load_spu = load.spu
        load._set_kW(kW)
        new_load_spu = load.spu
        bus._set_spu(old_load_spu, new_load_spu)

    def set_kvar(self, load_name: str, kvar: float):
        """
        sets a new kvar for the given Load.
        Updates Load.spu, Load.ppu, Load.qpu, and Bus.sum_spu
        """
        load = self.loads.get_ckt_element(load_name)
        bus = self.bus.get_ckt_element(load.bus_name)
        old_load_spu = load.spu
        load._set_kvar(kvar)
        new_load_spu = load.spu
        bus._set_spu(old_load_spu, new_load_spu)

    def get_lines_tx_idx_matrix(self):
        """
        n x 1 matrix of tx bus indices. Indexed by line index,
        which is the same value as in opendss
        """
        idx_matrix = np.zeros(self.lines.num_elements)
        for idx, line_name in self.lines._idx_to_name_dict.items():
            tx_bus = self.lines.get_element(line_name).tx
            idx_matrix[idx] = self.buses.get_idx(tx_bus)
        return idx_matrix

    def get_lines_rx_idx_matrix(self):
        """
        n x 1 matrix of rx bus indices. Indexed by line index,
        which is the same value as in opendss
        """
        idx_matrix = np.zeros(self.lines.num_elements)
        for idx, line_name in self.lines._idx_to_name_dict.items():
            rx_bus = self.lines.get_element(line_name).rx
            idx_matrix[idx] = self.buses.get_idx(rx_bus)
        return idx_matrix

    def get_spu_matrix(self) -> np.ndarray:
        """
        3 x n matrix of complex spu, columns indexed by bus index
        """
        spu_matrix = np.zeros((self.buses.num_elements, 3), dtype=complex)
        for load in self.loads.get_elements():
            bus_idx = self.buses.get_idx(load.related_bus)
            spu_matrix[bus_idx] += load.spu
        return spu_matrix.transpose()

    def get_cappu_matrix(self) -> np.ndarray:
        """
        3 x n matrix of complex spu, columns indexed by bus index
        """
        cappu_matrix = np.zeros((self.buses.num_elements, 3), dtype=complex)
        for cap in self.capacitors.get_elements():
            bus_idx = self.buses.get_idx(cap.related_bus)
            cappu_matrix[bus_idx] += cap.cappu
        return cappu_matrix.transpose()

    def get_aPQ_matrix(self) -> np.ndarray:
        """
        3 x n matrix of all load.aPQ_p, aggregated by phase on bus,
        columns indexed by bus
        """
        return self._get_zip_val_matrix('aPQ_p')

    def get_aI_matrix(self) -> np.ndarray:
        """
        3 x n matrix of all load.aPQ_p, aggregated by phase on bus,
        columns indexed by bus
        """
        return self._get_zip_val_matrix('aI_p')

    def get_aZ_matrix(self) -> np.ndarray:
        """
        3 x n matrix of all load.aPQ_p, aggregated by phase on bus,
        columns indexed by bus
        """
        return self._get_zip_val_matrix('aZ_p')

    def get_wpu_matrix(self) -> np.ndarray:
        """
        3 x n matrix of all wpu, columns indexed by bus
        Currently set to all zeros.
        TODO: Implement logic to set this as needed.
        """
        return np.zeros((3, self.buses.num_elements))

    def get_vvcpu_matrix(self) -> np.ndarray:
        """
        3 x n matrix of all wpu, columns indexed by bus
        Currently set to all zeros.
        TODO: Implement logic to set this as needed.
        """
        return np.zeros((3, self.buses.num_elements))

    def _get_zip_val_matrix(self, zip_param=str) -> np.ndarray:
        """
        3 x n matrix of all load.zip_param, aggregated by phase on bus,
        columns indexed by bus
        Zip_params can take any values in
        {'aPQ_p', 'aI_p', 'aZ_p','aPQ_q', 'aI_q', 'aZ_q'}
        """
        param_matrix = np.zeros((self.buses.num_elements, 3))
        for load in self.loads.get_elements():
            load_bus = load.related_bus
            load_ph_matrix = parse_phase_matrix(load.phases)
            bus_idx = self.buses.get_idx(load_bus)
            param_matrix[bus_idx] += load_ph_matrix * getattr(load, zip_param)
        return param_matrix.transpose()

    def _assign_to_buses(self, ckt_element_group):
        """
        For all elements in the ckt_element_group, gives the bus
        associated with CircuitElement.related_bus a pointer to the element
        """
        for ele in ckt_element_group.get_elements():
            bus = self.buses.get_element(ele.related_bus)
            element_list_ptr = f'{ele.__class__.__name__}s'.lower()
            if not getattr(bus, element_list_ptr):
                setattr(bus, element_list_ptr, [])
            getattr(bus, element_list_ptr).append(ele)
