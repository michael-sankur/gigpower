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
        self._init_KVL_matrices()
        self._init_KCL_matrices()
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

    def _init_KVL_matrices(self):
        """
        set self.b_KVL and self.G_KVL
        adapted from
        https://github.com/msankur/LinDist3Flow/blob/vectorized/20180601/PYTHON/lib/compute_SBKVL_matrices.py
        """
        tf_bus = self.circuit.transformers.get_bus_ph_matrix()
        tf_lines = self.circuit.transformers.get_num_lines_x_ph()
        vr_lines = self.circuit.voltage_regulators.get_num_lines_x_ph()
        nline = self.circuit.lines.num_elements
        nnode = self.circuit.buses.num_elements
        X_matrix = self.circuit.lines.get_X_matrix()
        R_matrix = self.circuit.lines.get_R_matrix()

        # ------- Residuals for KVL across line (m,n) ----------
        self.b_KVL = np.zeros((2*3*(nline) + 2*tf_lines + 2*2*vr_lines, 1))

        G_KVL = np.zeros((2*3*(nline) + 2*tf_lines + 2*2*vr_lines,
                         2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))

        for ph in range(3):
            for line in range(nline):  # line = line index
                line_ele = self.circuit.lines.get_element(line)
                bus1_idx = self.circuit.buses.get_idx(line_ele.tx)
                bus2_idx = self.circuit.buses.get_idx(line_ele.rx)
                bus1_phases = line_ele.phase_matrix
                if bus1_phases[ph] == 1:
                    G_KVL[2*ph*nline + 2*line][2*(nnode)*ph + 2*(bus1_idx)] = 1 #A_m
                    G_KVL[2*ph*nline + 2*line][2*(nnode)*ph + 2*(bus2_idx)] = -1 #A_n
                    G_KVL[2*ph*nline + 2*line+1][2*(nnode)*ph + 2*(bus1_idx) + 1] = 1 #B_m
                    G_KVL[2*ph*nline + 2*line+1][2*(nnode)*ph + 2*(bus2_idx) + 1] = -1 #B_n

                    G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 2*line] = -R_matrix[line][ph*3] * bus1_phases[0] #C_mn for a
                    G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 2*line + 1] = X_matrix[line][ph*3] * bus1_phases[0] #D_mn for a
                    G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 2*nline + 2*line] = -R_matrix[line][ph*3 + 1] * bus1_phases[1] #C_mn for b
                    G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 2*nline + 2*line + 1] = X_matrix[line][ph*3 + 1] * bus1_phases[1] #D_mn for b
                    G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 4*nline + 2*line] = -R_matrix[line][ph*3 + 2] * bus1_phases[2] #C_mn for c
                    G_KVL[2*ph*nline + 2*line][2*3*(nnode) + 4*nline + 2*line + 1] = X_matrix[line][ph*3 + 2] * bus1_phases[2] #D_mn for c

                    G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 2*line] = -X_matrix[line][ph*3] * bus1_phases[0] #C_mn for a
                    G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 2*line + 1] = -R_matrix[line][ph*3] * bus1_phases[0] #D_mn for a
                    G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 2*nline + 2*line] = -X_matrix[line][ph*3 + 1] * bus1_phases[1] #C_mn for b
                    G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 2*nline + 2*line + 1] = -R_matrix[line][ph*3 + 1] * bus1_phases[1] #D_mn for b
                    G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 4*nline + 2*line] = -X_matrix[line][ph*3 + 2] * bus1_phases[2] #C_mn for c
                    G_KVL[2*ph*nline + 2*line+1][2*3*(nnode) + 4*nline + 2*line + 1] = -R_matrix[line][ph*3 + 2] * bus1_phases[2] #D_mn for c
                else:
                    G_KVL[2*ph*nline + 2*line][2*(nnode)*3 + 2*ph*nline + 2*line] = 1 #C_mn
                    G_KVL[2*ph*nline + 2*line+1][2*(nnode)*3 + 2*ph*nline + 2*line+1] = 1 #D_mn

        #------- Residuals for Transformer KVL ----------

        line_idx_tf = range(0, tf_lines)
        kvl_count = 0
        for tfbs in range(len(tf_bus[0])):
            for ph in range(0, 3):
                if tf_bus[ph + 2, tfbs] != 0:
                    line = line_idx_tf[kvl_count]
                    G_KVL[2*3*nline + 2*line][2*nnode*ph + 2*int(tf_bus[0, tfbs])] = 1 #A_m
                    G_KVL[2*3*nline + 2*line][2*nnode*ph + 2*int(tf_bus[1, tfbs])] = -1 #A_n
                    G_KVL[2*3*nline + 2*line+1][2*nnode*ph + 2*int(tf_bus[0, tfbs]) + 1] = 1 #B_m
                    G_KVL[2*3*nline + 2*line+1][2*nnode*ph + 2*int(tf_bus[1, tfbs]) + 1] = -1 #B_n
                    kvl_count += 1

        self.G_KVL = G_KVL

    def _init_KCL_matrices(self):
        der, capacitance = 0, 0
        tf_bus = self.circuit.transformers.get_bus_ph_matrix()
        vr_bus = self.circuit.voltage_regulators.get_bus_ph_matrix()
        tf_lines = self.circuit.transformers.get_num_lines_x_ph()
        tf_no = self.circuit.transformers.num_elements
        vr_lines = self.circuit.voltage_regulators.get_num_lines_x_ph()
        vr_count = self.circuit.voltage_regulators.num_elements
        nnode = self.circuit.buses.num_elements
        nline = self.circuit.lines.num_elements

        # Line Indices Associated with Voltage Regulators and Transformers
        line_in_idx_vr = range(0, 2*vr_lines, 2)
        line_out_idx_vr = range(1, 2*vr_lines, 2)

        line_idx_tf = range(0, tf_lines)

        load_kw = self.circuit.loads.get_ppu_matrix()
        load_kvar = self.circuit.loads.get_qpu_matrix()
        caparr = self.circuit.capacitors.get_cappu_matrix()

        # ----------Residuals for KCL at a bus (m) ----------

        # Zip Parameters
        # Load
        beta_S = self.circuit.loads.aPQ_p
        beta_I = self.circuit.loads.aI_p
        beta_Z = self.circuit.loads.aZ_p

        # Capacitors
        gamma_S = self.circuit.capacitors.aPQ_p
        gamma_I = self.circuit.capacitors.aI_p
        gamma_Z = self.circuit.capacitors.aZ_p

        H = np.zeros((2*3*(nnode-1), 2*3*(nnode+nline) + 2*tf_lines + \
                      2*2*vr_lines, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
        g = np.zeros((2*3*(nnode-1), 1, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
        b = np.zeros((2*3*(nnode-1), 1, 1))

        # --------------- Quadratic Terms -----------------
        for ph in range(3):
            if ph == 0:  # set nominal voltage based on phase
                A0 = 1
                B0 = 0
            elif ph == 1:
                A0 = -1/2
                B0 = -1 * np.sqrt(3)/2
            elif ph == 2:
                A0 = -1/2
                B0 = np.sqrt(3)/2
            for k2 in range(1, nnode):  # skip slack bus
                bus = self.circuit.buses.get_element(k2)
                in_lines = self.circuit.lines.reverse_adj[bus.__name__]  # upstream buses
                out_lines = self.circuit.lines.adj[bus.__name__]  # downstream buses
                for cplx in range(2):
                    idxbs = k2
                    if cplx == 0:
                        load_val = load_kw[ph, idxbs]
                        cap_val = 0
                    else:
                        load_val = load_kvar[ph, idxbs]  # WPU HERE !!
                        cap_val = caparr[ph][idxbs]
                    # gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))]) # some derivatives
                    hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
                                            [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])
                    available_phases = bus.phase_matrix  # phase array at specific bus
                    if available_phases[ph] == 1:                 # quadratic terms
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2] = \
                                                -load_val * (beta_Z + (0.5 * beta_I * hessian_mag[0][0])) + \
                                                cap_val * (gamma_Z + (0.5 * gamma_I * hessian_mag[0][0]))  # TE replace assignment w/ -load_val * beta_Z; #a**2
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2 + 1] = \
                            -load_val * (beta_Z + (0.5 * beta_I * hessian_mag[1][1])) + \
                            cap_val * (gamma_Z + (0.5 * gamma_I * hessian_mag[1][1]))  # TE replace assignment w/ -load_val * beta_Z; #b**2
                            # H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1] = -load_val * beta_I * hessian_mag[0][1] / 2 #remove for TE
                            # H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2] =  -load_val * beta_I * hessian_mag[1][0] / 2 #remove for TE

                    for i in range(len(in_lines)):  # in lines
                        line_idx = self.circuit.lines.get_idx((in_lines[i], bus.__name__))
                        if available_phases[ph] == 1:
                            if cplx == 0:  # real residual
                                # A_m and C_lm
                                H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx] = 1/2
                                H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2] = 1/2
                                # B_m and D_lm
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = 1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = 1/2
                            if cplx == 1:  # imaginary residual
                                # #A_m, D_lm
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = -1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = -1/2
                                #B_m and C_lm
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx] = 1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = 1/2

                    for j in range(len(out_lines)):  # out lines
                        line_idx = self.circuit.lines.get_idx((bus.__name__, out_lines[j]))
                        if available_phases[ph] == 1:
                            if cplx == 0: #real residual
                                #A_m and C_mn
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx] = -1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2] = -1/2
                                #B_m and D_mn
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = -1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = -1/2
                            if cplx == 1: #imaginary residual
                                #A_m and D_mn
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1]= 1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = 1/2
                                #C_m and B_mn
                                H[2*ph*(nnode-1) + (k2-1)*2+cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx] = -1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = -1/2

        # ----------------------- Transformer KCL -----------------------
        count_tf = 0
        count_tf2 = 0
        for i in range(tf_no):
            for ph in range(0,3):
                k2 = int(tf_bus[1, i]) #in bus index of transformer [out bus: a0] --line-a0-a1-- [in bus: a1]
                if k2 == 0: #if source bus, need to skip line
                    count_tf += 1
                if k2 != 0 and tf_bus[ph + 2, i] != 0: #if not source bus, perform KCL
                    line_idx = line_idx_tf[count_tf]
                    #A_m and C_lm
                    H[2*ph*(nnode-1) + (k2-1)*2 ][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*line_idx] = 1/2
                    H[2*ph*(nnode-1) + (k2-1)*2 ][2*3*(nnode+nline) + 2*line_idx][2*(nnode)*ph + 2*k2] = 1/2
                    #B_m and D_lm
                    H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*line_idx + 1] = 1/2
                    H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = 1/2

                    #A_m, D_lm
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*line_idx + 1] = -1/2
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*3*(nnode+nline) + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = -1/2
                    #B_m and C_lm
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*line_idx] = 1/2
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*3*(nnode+nline) + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = 1/2
                    count_tf += 1 #go to next line

        for j in range(tf_no): #fill in H for the outlines
            for ph in range(0,3):
                k2 = int(tf_bus[0, j]) #out bus index of transformer
                if k2 == 0:
                    count_tf2 += 1
                if k2 != 0 and tf_bus[ph + 2, j] != 0:
                    line_idx = line_idx_tf[count_tf2]
                    #real residual
                    #A_m and C_mn
                    H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*line_idx] = -1/2
                    H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*line_idx][2*(nnode)*ph + 2*k2] = -1/2
                    #B_m and D_mn
                    H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*line_idx + 1] = -1/2
                    H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = -1/2

                    #imaginary residual
                    #A_m and D_mn
                    H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*line_idx + 1]= 1/2
                    H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*3*(nnode+nline) + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = 1/2
                    #C_m and B_mn
                    H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*line_idx] = -1/2
                    H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*3*(nnode+nline) + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = -1/2
                    count_tf2+=1

        # ----------------------- Voltage Regulator KCL -----------------------
        count_vr = 0
        count_vr2 = 0
        if vr_count > 0:
            for i in range(vr_count): #in lines
                for ph in range(0,3):
                    k2 = int(vr_bus[1, i])
                    if k2 == 0:
                        count_vr += 1
                    if k2 != 0 and vr_bus[ph + 2, i] != 0:
                        line_idx = line_in_idx_vr[count_vr]
                        #real residual
                        #A_m and C_lm
                        H[2*ph*(nnode-1) + (k2-1)*2 + 0][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx] = 1/2
                        H[2*ph*(nnode-1) + (k2-1)*2 + 0][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx][2*(nnode)*ph + 2*k2] = 1/2
                        #B_m and D_lm
                        H[2*ph*(nnode-1) + (k2-1)*2 + 0][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1] = 1/2
                        H[2*ph*(nnode-1) + (k2-1)*2 + 0][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = 1/2
                        #imaginary residual
                        # #A_m, D_lm
                        H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1] = -1/2
                        H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = -1/2
                        #B_m and C_lm
                        H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx] = 1/2
                        H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = 1/2
                        count_vr += 1

            for j in range(vr_count): #out lines
                for ph in range(0,3):
                    k2 = int(vr_bus[0, j])
                    if k2 == 0:
                        count_vr2 += 1
                    if k2 != 0 and vr_bus[ph + 2, j] != 0:
                        line_idx = line_out_idx_vr[count_vr2]
                        #real residual
                        #A_m and C_mn
                        H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx] = -1/2
                        H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx][2*(nnode)*ph + 2*k2] = -1/2
                        #B_m and D_mn
                        H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1] = -1/2
                        H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = -1/2
                        #imaginary residual
                        #A_m and D_mn
                        H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1] = 1/2
                        H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = 1/2
                        #C_m and B_mn
                        H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx] = -1/2
                        H[2*ph*(nnode-1) + (k2-1)*2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = -1/2
                        count_vr2 += 1


        # ----------------------- Linear Term & Constant Term -----------------------
        for ph in range(0,3):
            for k2 in range(1, nnode):
                for cplx in range(0,2):
                    bus = self.circuit.buses.get_element(k2)
                    available_phases = bus.phase_matrix #phase array at specific bus
                    idxbs = k2
                    if cplx == 0:
                        load_val = load_kw[ph][idxbs]
                    else:
                        load_val = load_kvar[ph][idxbs]

                    # Linear terms
                    g_temp = np.zeros(2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines)
                    if available_phases[ph] == 0: #if phase does not exist
                        g_temp[2*(ph)*nnode + 2*k2 + cplx] = 1
                    g[2*(nnode-1)*ph + 2*(k2-1) + cplx, 0,:] = g_temp

                    # Constant terms
                    if cplx == 0:
                        if der.real != 0:
                            b_factor = der.real
                            b_factor = 0
                        else:
                            b_factor = 0
                    elif cplx == 1:
                        if capacitance != 0 or der.imag != 0:
                            b_factor = capacitance - der.imag
                            b_factor = 0
                        else:
                            b_factor = caparr[ph][k2] # WPU HERE !!
                    else:
                        b_factor = 0

                    if available_phases[ph] == 0: #if phase does not exist
                        b_temp = 0
                    else:
                        b_temp = (-load_val * beta_S) + (b_factor * gamma_S)
                    b[2*(nnode-1)*ph + 2*(k2-1) + cplx][0][0] = b_temp

        self.H, self.g, self.b = H, g, b


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
