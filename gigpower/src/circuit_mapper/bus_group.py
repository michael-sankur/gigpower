from . circuit_element_group import CircuitElementGroup
from . bus import Bus


class BusGroup(CircuitElementGroup):
    dss_module_name = 'Bus'
    ele_class = Bus

    def _collect_names(self, dss):
        """ Override super() to use dss.Circuit.AllBusNames()"""
        self._names = dss.Circuit.AllBusNames()
        self._populate_name_idx_dicts()
