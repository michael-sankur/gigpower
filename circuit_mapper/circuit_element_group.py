from utils import parse_phase_matrix
import pandas as pd
import numpy as np
from typing import Tuple, Union
from circuit_element import CircuitElement


class CircuitElementGroup():
    def __init__(self, dss, **kwargs):
        self.num_elements = 0
        self._name_to_object_dict = {}
        self._collect_names(dss, **kwargs)
        self._collect_elements(dss, **kwargs)

    def _collect_elements(self, dss, **kwargs):
        dss_module = getattr(dss, f'{self.__class__.dss_module_name}')
        ele_class = self.__class__.ele_class

        # specify the dss method that sets the active element
        dss_set_active = dss_module.Name
        if self.__class__.__name__ == 'BusGroup':
            dss_set_active = dss.Circuit.SetActiveBus  # special case for Buses
        for name in self._names:
            dss_set_active(name)  # set as active element
            ele = ele_class(name, dss, **kwargs)  # create element
            self.add_element(ele)

    def _collect_names(self, dss, **kwargs):
        dss_module = getattr(dss, f'{self.__class__.dss_module_name}')
        self._names = dss_module.AllNames()  # preserve index order in opendss
        self._populate_name_idx_dicts()

    def _populate_name_idx_dicts(self):
        """ only call AFTER self._names is set! """
        self._name_to_idx_dict = {name: idx for idx, name in enumerate(self._names)}
        self._idx_to_name_dict = {idx: name for idx, name in enumerate(self._names)}

    def all_names(self):
        """ returns a View over all names in Group"""
        return self._name_to_idx_dict.keys()

    def get_idx(self, obj: Union[str, CircuitElement]) -> int:
        """
        return the index of the object within the Group
        param obj: the object's __name__, or the object itself
        """
        if not isinstance(obj, str):
            name = obj.__name__
        else:
            name = obj
        return int(self._name_to_idx_dict[name])

    def get_name(self, idx: int):
        """ return the name of the object given its Group idx"""
        return self._idx_to_name_dict[idx]

    def get_phase_matrix(self) -> np.ndarray:
        """
        3 x n phase matrix of 1's where phases are present, 0's otherwise
        columns indexed by element index, which is the same as in opendss
        """
        return self._get_attr_by_idx('phase_matrix', 'col')

    def get_element(self, key: Union[str, int, Tuple[str, str]]):
        """ Returns an element given a name, index, or tuple of (tx_name, rx_name)"""
        key_error_msg = "Invalid key. Key must be str, int, or (tx, rx) tuple."
        if isinstance(key, str):
            return self._name_to_object_dict[key]
        elif isinstance(key, int):
            return self._name_to_object_dict[self._idx_to_name_dict[key]]
        elif isinstance(key, tuple):
            try:
                return self._key_to_element_dict[key]()
            except KeyError:
                print(key_error_msg)
        else:
            raise KeyError(key_error_msg)

    def get_elements(self):
        """ returns an iterable View over all elements in the Group"""
        return self._name_to_object_dict.values()

    def add_element(self, ele):
        """
        adds a CircuitElement to this group
        with index = current value of self.num_elements, making
        this element the highest_index element in the group
        """
        if not isinstance(ele, self.__class__.ele_class):
            raise ValueError(f"Cannot add element {ele.__class__} to group {self.__class__}")
        if ele.__name__ not in self._names:
            self._names.append(ele.__name__)
            self._name_to_idx_dict[ele.__name__] = self.num_elements
            self._idx_to_name_dict[self.num_elements] = ele.__name__
        self._name_to_object_dict[ele.__name__] = ele
        self.num_elements = len(self._names)

    def _get_attr_by_idx(self, attr: str, orient='row') -> np.ndarray:
        """
        helper method to get an n x ? matrix of ele.attr values, indexed by
        the element index (same as self._names order, and the index order in
        opendss)
        param attr: name of the element param
        param orient: 'row' for n x ? matrix indexed by ele index, 'col' for its
        transpose
        """

        val = getattr(self.get_element(0), attr)

        try:
            attr_size = val.size  # np.ndarray.size
        except AttributeError:
            try:
                attr_size = len(val)  # len(List)
            except TypeError:
                attr_size = 1  # a scalar

        return_matrix = np.zeros((self.num_elements, attr_size))

        for ele, idx in self._name_to_idx_dict.items():
            obj = self.get_element(ele)
            return_matrix[idx] = getattr(obj, attr)

        if orient == 'col':
            return return_matrix.transpose()

        return return_matrix

    def _get_attr_by_bus(self, attr: str, orient='row') -> np.ndarray:
        """
        helper method to get a num_buses x ? matrix of ele.attr values, summed
        over Buses, indexed by Bus index
        param attr: name of the element param
        param orient: 'row' for n_buses x ? matrix indexed by bus index, 'col'
        for its transpose
        """
        val = getattr(self.get_element(0), attr)

        try:
            attr_size = val.size  # np.ndarray.size
        except AttributeError:
            try:
                attr_size = len(val)  # len(List)
            except TypeError:
                attr_size = 1  # a scalar

        return_matrix = np.zeros((self.buses.num_elements, attr_size))

        for ele, idx in self._name_to_idx_dict.items():
            obj = self.get_element(ele)
            bus_idx = self.buses.get_idx(obj.related_bus)
            return_matrix[bus_idx] += getattr(obj, attr)

        if orient == 'col':
            return return_matrix.transpose()

        return return_matrix
