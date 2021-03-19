# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create NR3 Solution class, a namespace for calculations used by nr3

from solution import Solution
import numpy as np


class SolutionNR3(Solution):

    # class variables set for all SolutionNR3 instances
    # TODO: If any of these need to be set by instance, move into self.__init__
    SLACKIDX = 0  # assume slack bus is at index 0
    VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])
    V0, I0 = None, None
    tolerance = 1e-9
    maxiter = 100

    def __init__(self, dss_fp: str):
        super().__init__(dss_fp)  # sets self.circuit
        self._init_XNR()
        self._init_slack_bus_matrices()
        # self._init_load_values()

    def _init_XNR(self):
        """
        adapted from
        https://github.com/msankur/LinDist3Flow/blob/vectorized/20180601/PYTHON/lib/basematrices.py
        """
        V0, I0 = self.__class__.V0, self.__class__.I0
        Vslack = self.__class__.VSLACK

        nline = self.circuit.lines.num_elements
        nnode = self.circuit.buses.num_elements
        tf_lines = self.circuit.transformers.get_num_lines_x_ph()
        vr_lines = self.circuit.voltage_regulators.get_num_lines_x_ph()

        # XNR order is bus voltages, line currents, transformer line currents,
        # voltage regulator line currents * 2
        XNR = np.zeros((2*3*(nnode + nline) + 2*tf_lines + 2*2*vr_lines, 1))

        # intialize node voltage portion of XNR
        if V0 is None or len(V0) == 0:
            for ph in range(3):
                for k1 in range(nnode):
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
            XNR[(2*3*nnode):] = 0.0*np.ones((6*nline + 2*tf_lines + 2*2*vr_lines, 1))

        # If initial I is given
        elif len(I0) != 0:
            for ph in range(3):
                for k1 in range(nline):
                    XNR[(2*3*nnode) + 2*ph*nline + 2*k1] = I0[ph, k1].real
                    XNR[(2*3*nnode) + 2*ph*nline + 2*k1+1] = I0[ph, k1].imag
            XNR[(2*3*nnode + 2*3*nline):] = np.zeros((len(XNR) - 2*3*nnode - 2*3*nline), 1)

        self.XNR = XNR

    def _init_slack_bus_matrices(self):
        """
        Initializes g_SB and b_SB
        adapted from
        https://github.com/msankur/LinDist3Flow/blob/vectorized/20180601/PYTHON/lib/basematrices.py
        """
        tf_lines = self.circuit.transformers.get_num_lines_x_ph()
        vr_lines = self.circuit.voltage_regulators.get_num_lines_x_ph()

        nline = self.circuit.lines.num_elements
        nnode = self.circuit.buses.num_elements

        Vslack = self.__class__.VSLACK

        # ------------ Slack Bus ------------------
        self.g_SB = np.zeros((6, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
        sb_idx = [0, 1, 2*nnode, (2*nnode)+1, 4*nnode, (4*nnode)+1]
        for i in range(len(sb_idx)):
            self.g_SB[i, sb_idx[i]] = 1

        self.b_SB = np.zeros((6, 1))
        for i in range(3):
            self.b_SB[2*i, 0] = Vslack[i].real
            self.b_SB[(2*i) + 1] = Vslack[i].imag

    def compute_SBKVL_matrices(self):
        """
        adapted from
        https://github.com/msankur/LinDist3Flow/blob/vectorized/20180601/PYTHON/lib/compute_SBKVL_matrices.py
        """
        tf_bus = self.circuit.transformers.get_bus_ph_matrix()
        vr_bus = self.circuit.voltage_regulators.get_bus_ph_matrix()
        tf_lines = self.circuit.transformers.get_num_lines_x_ph()
        vr_lines = self.circuit.voltage_regulators.get_num_lines_x_ph()
        tf_count = self.circuit.transformers.num_elements
        vr_no = self.circuit.voltage_regulators.num_elements
        gain = self.circuit.voltage_regulators.get_gain_matrix()

        nline = self.circuit.lines.num_elements
        nnode = self.circuit.buses.num_elements
        Sbase = self.circuit.Sbase

        # ---------- Resistance and Reactance Matrix for lines ---------------
        R_matrix = circuit.lines.get_R_matrix()
        X_matrix = circuit.lines.get_X_matrix()



    #     X = np.reshape(XNR, (2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines, 1))



    # #------- Residuals for KVL across line (m,n) ----------

    # G_KVL = np.zeros((2*3*(nline) + 2*tf_lines + 2*2*vr_lines,
    #                   2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))

    # for ph in range(0, 3):
    #     for line in range(len(dss.Lines.AllNames())):
    #         dss.Lines.Name(dss.Lines.AllNames()[line])  # set the line
    #         bus1 = dss.Lines.Bus1()
    #         bus2 = dss.Lines.Bus2()
    #         pattern = r"(\w+)\."
    #         try:
    #             bus1_idx = dss.Circuit.AllBusNames().index(
    #                 re.findall(pattern, bus1)[0])  # get the buses of the line
    #             bus2_idx = dss.Circuit.AllBusNames().index(
    #                 re.findall(pattern, bus2)[0])
    #         except:
    #             pattern = r"(\w+)"
    #             bus1_idx = dss.Circuit.AllBusNames().index(
    #                 re.findall(pattern, bus1)[0])  # get the buses of the line
    #             bus2_idx = dss.Circuit.AllBusNames().index(
    #                 re.findall(pattern, bus2)[0])
    #         # the buses on a line should have the same phase
    #         b1, _ = dss.CktElement.BusNames()
    #         if not dss.Lines.IsSwitch():
    #             # identifies which phase is associated with the bus (which is the same as the line)
    #             bus1_phases = identify_bus_phases(b1)
    #         else:
    #             # assume a switch is 3-phase
    #             bus1_phases = [1, 1, 1]
    #         if bus1_phases[ph] == 1:
    #             G_KVL[2*ph*nline + 2*line][2 *
    #                                        (nnode)*ph + 2*(bus1_idx)] = 1  # A_m
    #             G_KVL[2*ph*nline + 2*line][2 *
    #                                        (nnode)*ph + 2*(bus2_idx)] = -1  # A_n
    #             G_KVL[2*ph*nline + 2*line+1][2 *
    #                                          (nnode)*ph + 2*(bus1_idx) + 1] = 1  # B_m
    #             G_KVL[2*ph*nline + 2*line+1][2 *
    #                                          (nnode)*ph + 2*(bus2_idx) + 1] = -1  # B_n

    #             # C_mn for a
    #             G_KVL[2*ph*nline + 2*line][2*3 *
    #                                        (nnode) + 2*line] = -R_matrix[line][ph*3] * bus1_phases[0]
    #             # D_mn for a
    #             G_KVL[2*ph*nline + 2*line][2*3 *
    #                                        (nnode) + 2*line + 1] = X_matrix[line][ph*3] * bus1_phases[0]
    #             G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 2*nline + 2*line] = - \
    #                 R_matrix[line][ph*3 + 1] * bus1_phases[1]  # C_mn for b
    #             G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 2*nline + 2*line +
    #                                        1] = X_matrix[line][ph*3 + 1] * bus1_phases[1]  # D_mn for b
    #             G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 4*nline + 2*line] = - \
    #                 R_matrix[line][ph*3 + 2] * bus1_phases[2]  # C_mn for c
    #             G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 4*nline + 2*line +
    #                                        1] = X_matrix[line][ph*3 + 2] * bus1_phases[2]  # D_mn for c

    #             G_KVL[2*ph*nline + 2*line+1][2*3 *
    #                                          (nnode) + 2*line] = -X_matrix[line][ph*3] * bus1_phases[0]  # C_mn for a
    #             G_KVL[2*ph*nline + 2*line+1][2*3 *
    #                                          (nnode) + 2*line + 1] = -R_matrix[line][ph*3] * bus1_phases[0]  # D_mn for a
    #             G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 2*nline + 2*line] = - \
    #                 X_matrix[line][ph*3 + 1] * bus1_phases[1]  # C_mn for b
    #             G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 2*nline + 2*line + 1] = - \
    #                 R_matrix[line][ph*3 + 1] * bus1_phases[1]  # D_mn for b
    #             G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 4*nline + 2*line] = - \
    #                 X_matrix[line][ph*3 + 2] * bus1_phases[2]  # C_mn for c
    #             G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 4*nline + 2*line + 1] = - \
    #                 R_matrix[line][ph*3 + 2] * bus1_phases[2]  # D_mn for c
    #         else:
    #             G_KVL[2*ph*nline + 2*line][2 *
    #                                        (nnode)*3 + 2*ph*nline + 2*line] = 1  # C_mn
    #             G_KVL[2*ph*nline + 2*line+1][2 *
    #                                          (nnode)*3 + 2*ph*nline + 2*line+1] = 1  # D_mn

    # #------- Residuals for Transformer KVL ----------

    # line_idx_tf = range(0, tf_lines)
    # kvl_count = 0
    # for tfbs in range(len(tf_bus[0])):
    #     for ph in range(0, 3):
    #         if tf_bus[ph + 2, tfbs] != 0:
    #             line = line_idx_tf[kvl_count]
    #             G_KVL[2*3*nline + 2*line][2*nnode*ph +
    #                                       2*int(tf_bus[0, tfbs])] = 1  # A_m
    #             G_KVL[2*3*nline + 2*line][2*nnode*ph +
    #                                       2*int(tf_bus[1, tfbs])] = -1  # A_n
    #             G_KVL[2*3*nline + 2*line+1][2*nnode*ph +
    #                                         2*int(tf_bus[0, tfbs]) + 1] = 1  # B_m
    #             G_KVL[2*3*nline + 2*line+1][2*nnode*ph +
    #                                         2*int(tf_bus[1, tfbs]) + 1] = -1  # B_n
    #             kvl_count += 1

    # b_kvl = np.zeros((2*3*(nline) + 2*tf_lines + 2*2*vr_lines, 1))

    # # ---------- Voltage Regulator -----------

    # H_reg = np.zeros((2*vr_lines, 2*3*(nnode+nline) + 2*tf_lines +
    #                   2*2*vr_lines, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
    # G_reg = np.zeros(
    #     (2*vr_lines, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))

    # #  voltage ratio: V_bus2 - gamma V_bus1 = 0
    # line_in_idx = range(0, 2*vr_lines, 2)
    # line_out_idx = range(1, 2*vr_lines, 2)

    # vr_counter = 0
    # for m in range(vr_count):
    #     bus1_idx = vr_bus[0, m]
    #     bus2_idx = vr_bus[1, m]
    #     for ph in range(0, 3):
    #         if vr_bus[ph + 2, m] != 0:
    #             # voltage gain: gamma*A_in = A_out
    #             # gamma * B_in = B_out
    #             # Negative not shown below because inserted into gain term
    #             G_reg[2*vr_counter][2*nnode*ph + 2*bus1_idx] = gain[m]  # A_in
    #             G_reg[2*vr_counter][2*nnode*ph + 2*bus2_idx] = 1  # A_out
    #             G_reg[2*vr_counter + 1][2*nnode*ph +
    #                                     2*bus1_idx + 1] = gain[m]  # B_in
    #             G_reg[2*vr_counter + 1][2*nnode *
    #                                     ph + 2*bus2_idx + 1] = 1  # B_out

    #             #conservation of power: V_bus1 (I_bus1,out)* -  V_bus2 (I_bus2,in)* = 0
    #             # A_1 * C_out + B_1 * D_out  - (A_2 * C_in + B_2 * D_in)
    #             # j(B_1 * C_out - A_1 * D_out) - j(B_2 * C_in - A_2 * D_in)

    #             #A_1 C_out
    #             H_reg[2*vr_counter][2*nnode*ph + 2*bus1_idx][2*3 *
    #                                                          (nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter]] = 1
    #             H_reg[2*vr_counter][2*3*(nnode+nline) + 2*tf_lines + 2 *
    #                                 line_out_idx[vr_counter]][2*nnode*ph + 2*bus1_idx] = 1
    #             #B_1 D_out
    #             H_reg[2*vr_counter][2*nnode*ph + 2*bus1_idx + 1][2*3 *
    #                                                              (nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter] + 1] = 1
    #             H_reg[2*vr_counter][2*3*(nnode+nline) + 2*tf_lines + 2 *
    #                                 line_out_idx[vr_counter] + 1][2*nnode*ph + 2*bus1_idx + 1] = 1
    #             #A_2 C_in
    #             H_reg[2*vr_counter][2*nnode*ph + 2*bus2_idx][2*3 *
    #                                                          (nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter]] = -1
    #             H_reg[2*vr_counter][2*3*(nnode+nline) + 2*tf_lines + 2 *
    #                                 line_in_idx[vr_counter]][2*nnode*ph + 2*bus2_idx] = -1
    #             #B_2 D_in
    #             H_reg[2*vr_counter][2*nnode*ph + 2*bus2_idx + 1][2*3 *
    #                                                              (nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter] + 1] = -1
    #             H_reg[2*vr_counter][2*3*(nnode+nline) + 2*tf_lines + 2 *
    #                                 line_in_idx[vr_counter] + 1][2*nnode*ph + 2*bus2_idx + 1] = -1

    #             #B_1 * C_out
    #             H_reg[2*vr_counter + 1][2*nnode*ph + 2*bus1_idx + 1][2*3 *
    #                                                                  (nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter]] = 1
    #             H_reg[2*vr_counter + 1][2*3*(nnode+nline) + 2*tf_lines + 2 *
    #                                     line_out_idx[vr_counter]][2*nnode*ph + 2*bus1_idx + 1] = 1
    #             # A_1 * D_out
    #             H_reg[2*vr_counter + 1][2*nnode*ph + 2*bus1_idx][2*3 *
    #                                                              (nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter] + 1] = -1
    #             H_reg[2*vr_counter + 1][2*3*(nnode+nline) + 2*tf_lines + 2 *
    #                                     line_out_idx[vr_counter] + 1][2*nnode*ph + 2*bus1_idx] = -1
    #             #B_2 * C_in
    #             H_reg[2*vr_counter + 1][2*nnode*ph + 2*bus2_idx + 1][2*3 *
    #                                                                  (nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter]] = -1
    #             H_reg[2*vr_counter + 1][2*3*(nnode+nline) + 2*tf_lines + 2 *
    #                                     line_in_idx[vr_counter]][2*nnode*ph + 2*bus2_idx + 1] = -1
    #             # A_2 * D_in
    #             H_reg[2*vr_counter + 1][2*nnode*ph + 2*bus2_idx][2*3 *
    #                                                              (nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter] + 1] = 1
    #             H_reg[2*vr_counter + 1][2*3*(nnode+nline) + 2*tf_lines +
    #                                     2*line_in_idx[vr_counter] + 1][2*nnode*ph + 2*bus2_idx] = 1
    #             vr_counter += 1

    # return X, g_SB, b_SB, G_KVL, b_kvl, H_reg, G_reg


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
