# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Circuit class to mirror a dss Circuit object
# used by Solution objects to solve powerflow

import numpy as np

from . bus_group import BusGroup
from . line_group import LineGroup
from . load_group import LoadGroup
from . capacitor_group import CapacitorGroup
from . transformer_group import TransformerGroup
from . voltage_regulator_group import VoltageRegulatorGroup

from . utils import parse_phase_matrix, set_zip_values_dss


class Circuit():

    # zip values are set at the Circuit class level so that 
    # all Circuit instances use the same zip values
    ZIP_V = np.asarray([0.10, 0.05, 0.85, 0.10, 0.05, 0.85, 0.80])
    aPQ_p, aI_p, aZ_p = ZIP_V[0:3]
    aPQ_q, aI_q, aZ_q = ZIP_V[3:6]
    min_voltage_pu = ZIP_V[6]

    def __init__(self, dss, Sbase=10**6):
        """
        initialize Circuit from a dss object
        Note that the Solution class runs 'redirect' on the dss file
        The Circuit does not call opendss functions directly
        """
        set_zip_values_dss(dss, Circuit.ZIP_V)
        self.Sbase = Sbase
        self.buses = BusGroup(dss)
        self.lines = LineGroup(dss, bus_group=self.buses)
        self.loads = LoadGroup(dss, bus_group=self.buses)
        self.capacitors = CapacitorGroup(dss, bus_group=self.buses)
        self.transformers = TransformerGroup(dss, bus_group=self.buses, line_group=self.lines)
        self.voltage_regulators = VoltageRegulatorGroup(dss, line_group=self.lines)

        # # set zip values according to class vairables Circuit.ZIP_V 
        # self.loads._set_zip_values(Circuit.ZIP_V)
        # self.capacitors._set_zip_values(Circuit.ZIP_V)
        self.dss = dss

        # self._assign_to_buses(self.loads)
        # self._assign_to_buses(self.capacitors)
        # self._assign_to_buses(self.voltage_regulators)
        # self._assign_to_buses(self.transformers)

    @classmethod
    def set_zip_values(cls, zip_V):
        """
        sets zip values for the Circuit class
        param zip_V: List or nd.array with 7 values
        [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min voltage pu]
        Note that zip values are set only on the Circuit class.
        """
        cls.ZIP_V = np.asarray(zip_V)
        cls.aZ_p, cls.aI_p, cls.aPQ_p = cls.ZIP_V[0:3]
        cls.aZ_q, cls.aI_q, cls.aPQ_q = cls.ZIP_V[3:6]
        cls.min_voltage_pu = cls.ZIP_V[6]

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

    def get_tx_idx_matrix(self):
        """
        n x 1 matrix of tx bus indices, for all Lines
        Indexed as follows:
        [0, len(Lines) - 1]: Lines
        [len(Lines), len(Transformers)- 1]: Transformers
        [len(Transformers), len(VoltageRegulators)- 1]: VoltageRegulators

        """
        # nlines = self.lines.num_elements
        tx_buses = self.lines.get_bus_ids('tx')
        #     tx_buses += self.transformers.get_bus_ids('tx')
        #     tx_buses += self.voltage_regulators.get_bus_ids('tx')
        # except AttributeError:
        #     pass
        return np.asarray([self.buses.get_idx(bus) for bus in tx_buses])

    def get_rx_idx_matrix(self):
        """
        n x 1 matrix of rx bus indices. Indexed by line index,
        which is the same value as in opendss
        """
        # nlines = self.lines.num_elements
        rx_buses = self.lines.get_bus_ids('rx')# [0: nlines]
        # try:
        #     rx_buses += self.transformers.get_bus_ids('rx')
        #     rx_buses += self.voltage_regulators.get_bus_ids('rx')
        # except AttributeError:
        #     pass
        return np.asarray([self.buses.get_idx(bus) for bus in rx_buses])

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
    
    def get_total_lines(self):
        """ returns number of Lines and Synthetic Lines"""
        return len(self.lines.all_names())

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
            param_matrix[bus_idx] += load_ph_matrix * getattr(Circuit, zip_param)
        return param_matrix.transpose()

    def _assign_to_buses(self, ckt_element_group):
        """
        For all elements in the ckt_element_group, gives the bus
        associated with CircuitElement.related_bus a pointer to the element
        """
        for ele in ckt_element_group.get_elements():
            bus = self.buses.get_element(ele.related_bus)
            element_list_ptr = f'{ele.__class__.__name__}s'.lower()
            try:
                getattr(bus, element_list_ptr)
            except(AttributeError):
                setattr(bus, element_list_ptr, [])
            getattr(bus, element_list_ptr).append(ele)
