from circuit_element_group import CircuitElementGroup
from bus import Bus


class BusGroup(CircuitElementGroup):
    dss_module_name = 'Bus'
    ele_class = Bus

    def _collect_names(self, dss):
        """ Override super() to use dss.Circuit.AllBusNames()"""
        self._names = dss.Circuit.AllBusNames()
        self._name_to_idx_dict = {
            name: idx for idx, name in enumerate(self._names)}
        self._idx_to_name_dict = {
            idx: name for idx, name in enumerate(self._names)}
