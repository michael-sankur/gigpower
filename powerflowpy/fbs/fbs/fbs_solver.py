from . network import Network
from . solution import Solution
from . fbs import fbs
import numpy as np # type: ignore
from . utils import init_from_dss

class FBS:
    def __init__(self, dss_path:str) -> None:
        self.network = init_from_dss(dss_path)

    def solve(self) -> None:
        self.solution = fbs(self.network)

    def set_load_kw(self, load:str, kw: float) -> None:
        self.network.set_load_kw(load, kw)

    def set_load_kvar(self, load : str, kvar: float) -> None:
        self.network.set_load_kvar(load, kvar)

    def get_all_bus_voltages(self) -> None:
        all_voltages = np.absolute(np.array(list(self.solution.V.values())).flatten())
        return all_voltages[all_voltages != 0]
