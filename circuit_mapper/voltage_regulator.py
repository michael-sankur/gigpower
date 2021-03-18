from line import Line, SyntheticLine
from circuit_element import CircuitElement
import numpy as np
from utils import parse_dss_bus_name


class VoltageRegulator(CircuitElement):
    dss_module_name = 'RegControls'

    def __init__(self, name: str, dss, line_group):
        super().__init__(name, dss)
        self._set_lines(line_group)
        # set this regcontrol and its transformer as active,
        # in order to get the tap number
        dss.RegControls.Name(self.__name__)
        dss.Transformers.Name(self.transformer_name)
        self.tap = dss.RegControls.TapNumber()
        self.set_gamma(self.tap)
        self.Itx = np.zeros((3,), dtype='float')
        # complex current entering/leaving the regControl node
        self.Ireg = np.zeros((3,), dtype='float')

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
        self.related_bus = parse_dss_bus_name(regControl_bus)
        self.regControl_bus = self.related_bus  # alias for related bus
        self.tx = parse_dss_bus_name(tx)
        self.key = (self.tx, self.related_bus)

    def _set_phases(self, dss):
        """
        Use the default CircuitElement method to set phases by regControl Bus
        """
        CircuitElement._set_phases(self, dss)

    def _set_lines(self, line_group):
        """
        Each voltage regulator is modeled with two lines:
        1. upstream: txBus --> regControlBus
        2. downstream: regControlBus--> rxBus
        opendss already has a Line for the downstream line, so it should
        be present in Circuit.lines
        Creates a SyntheticLine for the upstream line, find the downstream line
        from the line_group,  and assign voltage regulators to both lines
        """
        self.upstream_line = SyntheticLine(self)
        self.downstream_line = self._find_downstream_line(line_group)
        self.upstream_line.add_voltage_regulator(self)
        self.downstream_line.add_voltage_regulator(self)

    def _find_downstream_line(self, line_group):
        for line in line_group.get_elements():
            if line.tx == self.regControl_bus:
                return line