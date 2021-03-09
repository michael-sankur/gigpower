from line import Line
import numpy as np
from utils import parse_dss_bus_name


class Transformer (Line):

    def __init__(self, dss):
        super().__init__(self, dss)
        # initialize this as a Line with 0 length, FZpu = zeroes(3x3)
        self.length = 0
        self.FZpu = np.zeros((3, 3), dtype='complex')
        # get kV, kVA
        dss.Transformers.Name(self.__name__)
        self.kV = dss.Transformers.kV()
        self.kVA = dss.Transformers.kVA()

    def _set_related_bus(self, dss):
        """
        Set tx and rx buses, set rated voltages based on tx and rx
        Set phases based on tx bus
        """
        dss.Transformers.Name(self.__name__)
        tx, rx = (parse_dss_bus_name(b) for b in dss.CktElement.BusNames())
        self.key = (tx, rx)
        upstream = self._get_rated_voltages(dss, 1)
        downstream = self._get_rated_voltages(dss, 2)
        self.rated_voltages = (upstream, downstream)
        self._set_phases_from_bus(tx)

    def _get_rated_voltages(self, dss, winding):
        """ --winding-- Set to 1 for upstream side, 2 for downstream side"""
        if dss.CktElement.NumPhases() == 1:  # 1 phase -> LN voltage
            return dss.Transformers.kV()
        return dss.Transformers.kV() / (3 ** 0.5)  # 2 or 3 phases -> LL voltage
