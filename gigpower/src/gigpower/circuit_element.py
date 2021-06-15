# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021

from . utils import parse_phase_matrix, parse_dss_bus_name
from typing import Union
from importlib import import_module
import numpy as np


class CircuitElement():
    def __init__(self, name, dss, **kwargs):
        self.__name__ = name
        self.phases = ''
        self.Sbase = 10**6
        self._set_related_bus(dss, **kwargs)
        self._set_base_vals(dss, **kwargs)
        self._set_phases(dss, **kwargs)

    def __str__(self) -> str:
        return f"{self.__class__}, {self.name}, {self.phases}"

    def get_ph_idx_matrix(self) -> np.ndarray:
        """
        return indices of active phases on element
        ex: self.phases = ['1', '2'] returns [0, 1]
        """
        return np.asarray([int(ph) - 1 for ph in self.phases])

    def _set_related_bus(self, dss):
        """
        Get the bus asssociated from this element.
        Set base vals, phases, and self.bus_name based on this bus
        """
        dss_module = import_module(f'opendssdirect.{self.__class__.dss_module_name}')
        dss_module.Name(self.__name__)  # set this element as active
        bus_name = dss.CktElement.BusNames()[0]  # this is usually the related bus
        self.related_bus = parse_dss_bus_name(bus_name)

    def _set_base_vals(self, dss, **kwargs):
        """ set Vbase, Ibase, and Zbase based self.related_bus"""
        dss.Circuit.SetActiveBus(self.related_bus)
        self.Vbase = dss.Bus.kVBase() * 1000
        self.Ibase = self.Sbase/self.Vbase
        self.Zbase = self.Vbase/self.Ibase

    def _set_phases(self, dss):
        """ set element's phases based on self.related_bus"""
        dss.Circuit.SetActiveBus(self.related_bus)
        self.phases = dss.Bus.Nodes()
        self.phase_matrix = parse_phase_matrix(self.phases)

    def _set_attr_val_by_phase(self, attr: str, value: Union[float, complex]):
        """
        sets self.attr to a 1x3 matrix, where phases present on self
        are set to value, 0 otherwise
        EX: self.phases = ['1', '3']
        self._set_attr_val_by_phase('ppu', .38) -> self.ppu = [.38, 0, .38]
        """
        setattr(self, attr, self.phase_matrix * value)

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

