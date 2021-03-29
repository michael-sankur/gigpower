from . line import Line, SyntheticLine
from . circuit_element import CircuitElement
import numpy as np
from . utils import parse_dss_bus_name, parse_dss_phases, parse_phase_matrix
from typing import Tuple

class VoltageRegulator(CircuitElement):
    dss_module_name = 'RegControls'

    def __init__(self, name: str, dss):
        super().__init__(name, dss)
        self._set_lines()
        # set this regcontrol and its transformer as active,
        # in order to get the tap number
        dss.RegControls.Name(self.__name__)
        dss.Transformers.Name(self.transformer_name)
        self.tap = dss.RegControls.TapNumber()
        self.set_gamma(self.tap)
        self.Itx = np.zeros((3,), dtype=complex)
        # complex current entering/leaving the regControl node
        self.Ireg = np.zeros((3,), dtype=complex)

    def set_gamma(self, tapNumber: int) -> None:
        """
        Takes a Tap Number from opendss and maps it to a voltage ratio.
        sets self.gamma to the result.
        There are 33 default integer Tap Numbers, from -16 to +16,
        corresponding to the defaults of 0.90 to 1.10 voltage ratio
        """
        VRMIN = .90
        VRMAX = 1.10
        # move range [-16,16] to [0,1]
        result = (tapNumber + 16) / 32
        # compress to VRMAX - VRMIN
        result *= VRMAX - VRMIN
        # move to [VRMIN, VRMAX]
        result += VRMIN
        # assign to self
        self.gamma = result

    def _set_related_bus(self, dss):
        """
        Overrides superclass method to get both tx and regControl Bus.
        Set base vals, phases, and self.bus_name based on this bus
        """
        dss.RegControls.Name(self.__name__)
        self.transformer_name = dss.RegControls.Transformer()
        dss.Transformers.Name(self.transformer_name)
        tx, regControl_bus = dss.CktElement.BusNames()  # get upstream, regcontrol buses
        self.tx, self.rx = parse_dss_bus_name(
            tx),  parse_dss_bus_name(regControl_bus)
        self.regControl_bus = self.rx
        self.related_bus = self.regControl_bus  # alias for related bus
        self.key = (self.tx, self.rx)

    def _set_phases(self, dss):
        """
        override super() to set name from CktElement.BusNames
        """
        dss.RegControls.Name(self.__name__)
        self.phases = parse_dss_phases(dss.CktElement.BusNames()[0])
        self.phase_matrix = parse_phase_matrix(self.phases)

    def _set_lines(self):
        """
        Each voltage regulator is modeled with two lines:
        1. txBus --> regControlBus
        2. regControlBus--> txBus
        Finds or creates this pair of SyntheticLines
        from the line_group,  and assigns voltage regulators to both lines
        Add the upstream line to the topology for the main line_group
        and to the key_to_element_dict
        """
        self.line_tx_to_reg = SyntheticLine(unique_key=False,
                                            name=self.__name__ + '_to_reg',
                                            key=self.key)
        self.line_reg_to_tx = SyntheticLine(unique_key=False,
                                            name=self.__name__ + '_to_tx', 
                                            key=self.key[-1::-1])
        # self.downstream_line = self._find_downstream_line(line_group)
        self.line_tx_to_reg.add_voltage_regulator(self)
        self.line_reg_to_tx.add_voltage_regulator(self)
        # self.downstream_line.add_voltage_regulator(self)

    # def _find_downstream_line(self, line_group):
    #     for line in line_group.get_elements():
    #         if line.tx == self.regControl_bus:
    #             return line
    #     return None
