# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Circuit class to mirror a dss Circuit object
# used by Solution objects to solve powerflow

from collections import OrderedDict
import numpy as np
import pandas as pd

from bus_group import BusGroup
from capacitor_group import CapacitorGroup
from line_group import LineGroup
from load_group import LoadGroup
from transformer_group import TransformerGroup
from voltage_regulator_group import VoltageRegulatorGroup


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

    def _assign_to_buses(self, ckt_element_group):
        """
        For all elements in the ckt_element_group, gives the bus
        associated with CircuitElement.bus_name a pointer to the element
        """
        for ele in ckt_element_group.elements():
            bus = self.buses.get_ckt_element(ele.bus_name)
            element_list_ptr = f'{ele.__class__.__name__}s'.lower()
            if not getattr(bus, element_list_ptr):
                setattr(bus, element_list_ptr, [])
            getattr(bus, element_list_ptr).append(ele)
