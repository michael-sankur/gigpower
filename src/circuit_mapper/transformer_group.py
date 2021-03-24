from . line_group import LineGroup
from . transformer import Transformer
import numpy as np


class TransformerGroup(LineGroup):
    dss_module_name, ele_class = 'Transformers', Transformer

    def _collect_names(self, dss, **kwargs):
        """ Override LineGroup._collect_names to exclude voltage regulators"""
        transformers = dss.Transformers.AllNames()
        vregs = dss.RegControls.AllNames()
        self._names = [n for n in transformers if n not in vregs]
        super()._populate_name_idx_dicts()
