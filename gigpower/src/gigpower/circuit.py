# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Circuit class to mirror a dss Circuit object
# used by Solution objects to solve powerflow

import numpy as np
import pandas as pd

from . bus_group import BusGroup
from . line_group import LineGroup
from . load_group import LoadGroup
from . capacitor_group import CapacitorGroup
from . transformer_group import TransformerGroup
from . voltage_regulator_group import VoltageRegulatorGroup

from . utils import parse_phase_matrix, set_zip_values_dss


class Circuit():
    @classmethod
    def _set_zip_values(cls, zip_V):
        """
        sets zip values for the Circuit class
        same method as Solution.set_zip_values, just private
        param zip_V: List or nd.array with 7 values
        [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min voltage pu]
        Note that zip values are set both on the Solution class and Circuit
        class. Users should set zip values via the Solution class. 
        """
        cls.ZIP_V = np.asarray(zip_V)
        cls.aZ_p, cls.aI_p, cls.aPQ_p = cls.ZIP_V[0:3]
        cls.aZ_q, cls.aI_q, cls.aPQ_q = cls.ZIP_V[3:6]
        cls.min_voltage_pu = cls.ZIP_V[6]

    def __init__(self, dss, Sbase=10**6):
        """
        initialize Circuit from a dss object
        Note that the Solution class runs 'redirect' on the dss file
        The Circuit does not call opendss functions directly
        """
        set_zip_values_dss(dss, Circuit.ZIP_V)
        self.Sbase = Sbase
        #: The Circuit's BusGroup
        self.buses = BusGroup(dss)
        #: The Circuit's LineGroup
        self.lines = LineGroup(dss, bus_group=self.buses)
        #: The Circuit's LoadGroup
        self.loads = LoadGroup(dss, bus_group=self.buses, zip_v=Circuit.ZIP_V)
        #: The Circuit's CapacitorGroup
        self.capacitors = CapacitorGroup(dss, bus_group=self.buses)
        #: The Circuit's TransformerGroup
        self.transformers = TransformerGroup(dss, bus_group=self.buses)
        #: The Circuit's VoltageRegulatorGroup
        self.voltage_regulators = VoltageRegulatorGroup(dss, bus_group=self.buses)

        # the main line group needs to be aware of transformers and voltage 
        # regulators. It can be queried for transformer and voltage regulator
        # indices and topology
        self.lines.transformers = self.transformers
        self.lines.voltage_regulators = self.voltage_regulators
        #: A pointer to the OpenDSS object corresponding to this Circuit
        self.dss = dss
        self._orient = 'rows'  # may be overwritten by Solution

    def set_kW(self, load_name: str, kW: float):
        """
        sets a new kW for the given Load.
        Updates Load.spu, Load.ppu, Load.qpu, and Bus.sum_spu
        """
        load = self.loads.get_element(load_name)
        bus = self.buses.get_element(load.related_bus)
        old_load_spu = load.spu
        load._set_kW(kW)
        new_load_spu = load.spu
        bus._set_spu(old_load_spu, new_load_spu)

    def set_kvar(self, load_name: str, kvar: float):
        """
        sets a new kvar for the given Load.
        Updates Load.spu, Load.ppu, Load.qpu, and Bus.sum_spu
        """
        load = self.loads.get_element(load_name)
        bus = self.buses.get_element(load.related_bus)
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
        tx_buses = self.lines.get_bus_ids('tx')
        try:
            tx_buses += self.transformers.get_bus_ids('tx')
            tx_buses += self.voltage_regulators.get_bus_ids('tx')
        except AttributeError:
            pass
        return np.asarray([self.buses.get_idx(bus) for bus in tx_buses])

    def get_rx_idx_matrix(self):
        """
        n x 1 matrix of rx bus indices. Indexed by line index,
        which is the same value as in opendss
        """
        rx_buses = self.lines.get_bus_ids('rx')
        try:
            rx_buses += self.transformers.get_bus_ids('rx')
            rx_buses += self.voltage_regulators.get_bus_ids('rx')
        except AttributeError:
            pass
        return np.asarray([self.buses.get_idx(bus) for bus in rx_buses])

    def _orient_switch(self, matrix):
        if self._orient == 'rows':
            return matrix
        elif self._orient == 'cols':
            return matrix.transpose()

    def get_spu_matrix(self) -> np.ndarray:
        """
        3 x n or n x 3 matrix of complex spu indexed by bus index
        """
        spu_matrix = self.loads.get_spu_matrix()
        return self._orient_switch(spu_matrix)

    def get_cappu_matrix(self) -> np.ndarray:
        """
        3 x n or n x 3 matrix of real cappu, columns indexed by bus index
        """
        cappu_matrix = self.capacitors.get_cappu_matrix()
        return self._orient_switch(cappu_matrix)

    def get_aPQ_matrix(self) -> np.ndarray:
        """
        3 x n or n x 3 matrix of all load.aPQ_p, aggregated by phase on bus,
        columns indexed by bus
        """
        matrix = self.loads._get_zip_val_matrix('aPQ_p')
        return self._orient_switch(matrix)

    def get_aI_matrix(self) -> np.ndarray:
        """
        3 x n or n x 3matrix of all load.aPQ_p, aggregated by phase on bus,
        columns indexed by bus
        """
        matrix = self.loads._get_zip_val_matrix('aI_p')
        return self._orient_switch(matrix)

    def get_aZ_matrix(self) -> np.ndarray:
        """
        3 x n or n x 3matrix of all load.aPQ_p, aggregated by phase on bus,
        columns indexed by bus
        """
        matrix = self.loads._get_zip_val_matrix('aZ_p')
        return self._orient_switch(matrix)

    def get_wpu_matrix(self) -> np.ndarray:
        """
        3 x n or n x 3matrix of all real wpu, columns indexed by bus
        Currently set to all zeros.
        TODO: Implement logic to set this as needed.
        """
        return self._orient_switch(np.zeros((self.buses.num_elements, 3),
                                   dtype=float))

    def get_total_lines(self):
        """ returns number of Lines transformers, and voltage regulators * 2"""
        total = self.lines.num_elements
        try:
            total += self.transformers.num_elements
        except AttributeError:
            pass
        try:
            total += self.voltage_regulators.num_elements
        except AttributeError:
            pass
        return total

    def get_nominal_bus_powers(self) -> pd.DataFrame:
        """ 3 x n or n x 3 matrix of total nominal powers by bus"""
        data = self.get_spu_matrix() - 1j * self.get_cappu_matrix()
        data_type = complex
        index = self.buses.all_names()
        cols = ['A', 'B', 'C']
        if self._orient == 'cols':
            # force a deep copy swap to avoid pointer issues
            temp = [_ for _ in cols]
            cols = [_ for _ in index]
            index = temp
        return pd.DataFrame(data=data, index=index, columns=cols, dtype=data_type)

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
