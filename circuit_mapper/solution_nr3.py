# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create NR3 Solution class, a namespace for calculations used by nr3

from solution import Solution
import numpy as np


class SolutionNR3(Solution):

    # class variables set for all SolutionNR3 instances
    # TODO: If any of these need to be set by instance, move into self.__init__
    SLACKIDX = 0
    VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])
    V0, I0 = None
    tolerance = 1e-9
    maxiter = 100

    def __init__(self, dss_fp: str):
        super().__init__(dss_fp)  # sets self.circuit
        self._init_XNR()
        self._init_load_values()

    def _init_XNR(self):
        """
        adapted from
        https://github.com/msankur/LinDist3Flow/blob/vectorized/20180601/PYTHON/lib/NR3.py
        """
        V0, I0 = self.__class__.V0, self.__class__.I0
        Vslack = self.__class__.VSLACK
        nnode = self.circuit.buses.num_elements
        nline = self.circuit.lines.num_elements
        XNR = np.zeros((2*3*(nnode + nline), 1))

        # intialize node voltage portion of XNR
        if V0 is None or len(V0) == 0:
            for ph in range(0, 3):
                for k1 in range(0, nnode):
                    XNR[2*ph*nnode + 2*k1] = Vslack[ph].real
                    XNR[2*ph*nnode + 2*k1+1] = Vslack[ph].imag

        # If initial V is given (usually from CVX)
        elif len(V0) != 0:
            for ph in range(3):
                for k1 in range(nnode):
                    XNR[2*ph*nnode + 2*k1] = V0[ph, k1].real
                    XNR[2*ph*nnode + 2*k1+1] = V0[ph, k1].imag

        # intialize line current portion of XNR
        if I0 is None or len(I0) == 0:
            for k1 in range(nnode):
                XNR[(2*3*nnode):] = 0.0*np.ones((6*nline, 1))

        elif len(I0) != 0:
            for ph in range(3):
                for k1 in range(nline):
                    XNR[(2*3*nnode) + 2*ph*nline + 2*k1] = I0[ph, k1].real
                    XNR[(2*3*nnode) + 2*ph*nline + 2*k1+1] = I0[ph, k1].imag

        self.XNR = XNR

    def _init_load_values(self):
        load_order_list = self._init_load_order_f()
        bus_load = np.zeros((3, self.nnode, 2))
        load_ph_arr = np.zeros((self.nload, 3))

        load_ph_arr_origin = np.zeros((self.nnode, max(load_order_list.values()), 3))
        bus_load_divide = np.zeros((3, self.nnode, 2))

        for load in range(self.nload):
            dss.Loads.Name(self.all_load_names[load])
            pattern =  r"(\w+)\."
            load_bus = re.findall(pattern, dss.CktElement.BusNames()[0])
            load_ph_arr_temp = [0, 0, 0]
            for i in range(1, 4):
                pattern = r"\.%s" % (str(i))
                load_ph = re.findall(pattern, dss.CktElement.BusNames()[0])
                if load_ph:
                    load_ph_arr_temp[i - 1] = 1
                    load_ph_arr[load, i - 1] = 1
            for j in range(max(load_order_list.values())):
                idxbs = dss.Circuit.AllBusNames().index(load_bus[0])
                if np.all(load_ph_arr_origin[idxbs, j,:] == [0, 0, 0]):
                    load_ph_arr_origin[idxbs, j, :] = load_ph_arr_temp
                    for i in range(len(load_ph_arr_temp)):
                        if load_ph_arr_temp[i] == 1:
                            #bus_load[i, idxbs, 0] += dss.Loads.kW() *1e3*1 / self.Sbase / sum(load_ph_arr_temp)
                            #bus_load[i, idxbs, 1] += dss.Loads.kvar()*1e3*1 / self.Sbase  / sum(load_ph_arr_temp)
                            bus_load_divide[i, idxbs, 0] = 1e3 / self.Sbase / sum(load_ph_arr_temp)
                            bus_load_divide[i, idxbs, 1] = 1e3 / self.Sbase  / sum(load_ph_arr_temp)
                    break
        return bus_load, bus_load_divide, load_ph_arr
