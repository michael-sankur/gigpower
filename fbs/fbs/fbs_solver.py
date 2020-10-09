from . network import Network
from . solution import Solution
from . fbs import fbs
import numpy as np
from . utils import init_from_dss

class FBS:
    def __init__(self, dss_path):
        self.network = init_from_dss(dss_path)

    def solve(self):
        self.solution = fbs(self.network)

    def set_load_kw(self, load, kw):
        self.network.set_load_kw(load, kw)

    def set_load_kvar(self, load, kvar):
        self.network.set_load_kvar(load, kvar)

    def get_all_bus_voltages(self):
        self.solution.calc_sV()