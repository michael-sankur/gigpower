# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021

from utils import parse_phases
from typing import List

class CircuitElement():
    def __init__(self, name, dss):
        self.__name__ = name
        self.phases = ''
        self.Sbase = 10**6
        self._set_related_bus(dss)

    def __str__(self) -> str:
        return f"{self.__class__}, {self.name}, {self.phases}"

    def get_phase_matrix(self) -> List[bool]:
        return parse_phases(self.phases)

    def _set_related_bus(self, dss):
        """
        Get the bus asssociated from this element.
        Set base vals, phases, and self.bus_name based on this bus
        """
        dss_module = __import__(f'opendssdirect.{self.__class__.dss_module_name}')
        dss_module.Name(self.__name__)  # set this element as active
        bus_name = dss.CktElement.BusNames()[0]  # this is usually the related bus
        self._set_base_vals(bus_name)
        self._set_phases_from_bus(bus_name)
        self.bus_name = bus_name.split('.')[0]

    def _set_base_vals(self, dss, bus_name: str):
        """ set Vbase, Ibase, and Zbase based on a certain bus"""
        dss.Circuit.SetActiveBus(bus_name)
        self.Vbase = dss.Bus.kVBase() * 1000
        self.Ibase = self.Sbase/self.Vbase
        self.Zbase = self.Vbase/self.Ibase

    def _set_phases_from_bus(self, dss, bus_name: str):
        """ set element's phases based on a certain bus"""
        # Handle the typical case where bus names tell you the phases, e.g. 'BusName.1.2.'
        if "." in bus_name:
            bus_name, *phases = bus_name.split('.')
        if not phases:  # if no phases are present in name, assume all 3 phases
            phases = ['1', '2', '3']
        self.phases = phases

    def get_spu_matrix(self, dss):
        pass

    def get_aPQ_matrix(self, dss):
        pass

    def get_aZ_matrix(self, dss):
        pass

    def get_aI_matrix(self, dss):
        pass

    def get_cappu_matrix(self, dss):
        pass

    def get_wpu_matrix(self, dss):
        pass

    def get_vvcpu_matrix(self, dss):
        pass
