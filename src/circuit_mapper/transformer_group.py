from . line_group import LineGroup
from . transformer import Transformer
import numpy as np


class TransformerGroup(LineGroup):
    dss_module_name, ele_class = 'Transformers', Transformer

    def _collect_names(self, dss, **kwargs):
        """
        Override LineGroup._collect_names to exclude voltage regulators
        Note: the dss config file must give voltage regulators the same name
        as their transformers, otherwise this parser won't work
        TODO: document this requirement
        """
        transformers = dss.Transformers.AllNames()
        vregs = dss.RegControls.AllNames()
        self._names = [n for n in transformers if n not in vregs]
        super()._populate_name_idx_dicts()
