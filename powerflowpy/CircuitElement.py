# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021

import numpy as np  # type: ignore
import pandas as pd  # type: ignore


class CircuitElement():
    def __init__(self, name):
        self.__name__ = name
        self.phases = ''

    def __str__(self) -> str:
        return f"{self.__class__}, {self.name}, {self.phases}"
