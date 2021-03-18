from utils import parse_phase_matrix
import pandas as pd
import numpy as np
from typing import Tuple, Union


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

    def get_idx(self, name: str) -> int:
        """ return the index of the object within the Group given its name"""
        return int(self._name_to_idx_dict[name])

    def get_name(self, idx: int):
        """ return the name of the object given its Group idx"""
        return self._idx_to_name_dict[idx]

    def get_phase_matrix(self) -> np.ndarray:
        """
        3 x n phase matrix of 1's where phases are present, 0's otherwise
        columns indexed by element index, which is the same as in opendss
        """
        phase_matrix = np.zeros((len(self._names), 3), dtype=int)
        for ele, idx in self._name_to_idx_dict.items():
            obj = self.get_element(ele)
            phase_matrix[idx] = parse_phase_matrix(obj.phases)
        return phase_matrix.transpose()

    def get_phase_df(self) -> pd.DataFrame:
        """ n x 3 dataframe indexed by object names """
        phase_matrix = self.get_phase_matrix()
        return pd.DataFrame(data=phase_matrix, index=self._names, columns=['A', 'B', 'C'])

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
