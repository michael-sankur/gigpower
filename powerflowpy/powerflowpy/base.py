from typing import Dict, List

class Network:
    """
    Base class that defines a network. 
    """
    def __init__(self, dss_file = None) -> None:
        self.buses = [] # an array of Busses in sorted order of increasing distance from source bus.
        self.lines = [] # an array of lines

class Bus:
    pass

class Line:
    pass
