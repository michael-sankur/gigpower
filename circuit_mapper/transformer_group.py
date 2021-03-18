from line_group import LineGroup
from transformer import Transformer
import numpy as np


class TransformerGroup(LineGroup):
    dss_module_name, ele_class = 'Transformers', Transformer

    def __init__(self, dss, bus_group):
        super().__init__(dss)
        self.buses = bus_group

    def _collect_names(self, dss):
        """ Override LineGroup._collect_names to exclude voltage regulators"""
        transformers = dss.Transformers.AllNames()
        vregs = dss.RegControls.AllNames()
        self._names = [n for n in transformers if n not in vregs]
        super()._populate_name_idx_dicts()

    def get_tf_bus_ph_matrix(self):
        """
        5 x num_transformers matrix, columns indexed by transformer index
        row 0: tx_bus_index
        row 1: rx_bus_index
        row 2 - 4: A, B, C phases
        bus_group: the BusGroup of the Circuit
        """
        tf_matrix = np.zeros((self.num_elements, 5), dtype=int)
        for tf in self.get_elements():
            tf_idx = self.get_idx(tf.__name__)
            tx_bus_name, rx_bus_name = tf.key
            tf_matrix[tf_idx][0:2] = np.asarray([self.buses.get_idx(
                tx_bus_name), self.buses.get_idx(rx_bus_name)])
            tf_matrix[tf_idx][2:] = tf.get_phase_matrix()
        return tf_matrix.transpose()
