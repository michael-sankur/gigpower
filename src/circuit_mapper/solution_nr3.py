# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create NR3 Solution class, a namespace for calculations used by nr3

from . solution import Solution
from . circuit import Circuit
import numpy as np


class SolutionNR3(Solution):
    # TODO: set a class-level NR3 tolerance value

    def __init__(self, dss_fp: str, **kwargs):
        super().__init__(dss_fp, **kwargs)  # sets self.circuit
        self._set_orient('cols')
        self._init_XNR()
        self._init_slack_bus_matrices()
        self._init_KVL_matrices()
        self._init_KCL_matrices()
        self._init_KVL_matrices_vregs()

    def _init_XNR(self):
        """
        adapted from
        https://github.com/msankur/LinDist3Flow/blob/vectorized/20180601/PYTHON/lib/basematrices.py
        written by @kathleenchang
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
            XNR[(2*3*nnode):] = 0.0*np.ones((6*nline + 2*tf_lines + 2*2* vr_lines, 1))

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
        written by @kathleenchang
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

    def _init_KVL_matrices(self):
        """
        set self.b_KVL and self.G_KVL
        copied from
        https://github.com/msankur/LinDist3Flow/blob/vectorized/20180601/PYTHON/lib/compute_SBKVL_matrices.py
        written by @kathleenchang
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

    def _init_KCL_matrices(self, der=0, capacitance=0):
        """
        set H, g, b
        copied from 20180601/PYTHON/lib/compute_KCL_matrices.py
        written by @kathleenchang
        """
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
        bp = self.circuit.buses.get_phase_matrix_dict()
        dss = self.dss

        # Zip Parameters
        # Load
        beta_S = Circuit.aPQ_p
        beta_I = Circuit.aI_p
        beta_Z = Circuit.aZ_p

        # Capacitors
        gamma_S = Circuit.aPQ_p
        gamma_I = Circuit.aI_p
        gamma_Z = Circuit.aZ_p

        H = np.zeros((2*3*(nnode-1), 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
        g = np.zeros((2*3*(nnode-1), 1, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
        b = np.zeros((2*3*(nnode-1), 1, 1))

        # --------------- Quadratic Terms -----------------
        for ph in range(0,3):
            if ph == 0: #set nominal voltage based on phase
                A0 = 1
                B0 = 0
            elif ph == 1:
                A0 = -1/2
                B0 = -1 * np.sqrt(3)/2
            elif ph == 2:
                A0 = -1/2
                B0 = np.sqrt(3)/2
            for k2 in range(1, len(dss.Circuit.AllBusNames())): #skip slack bus
                dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k2]) #set the bus
                bus_name = dss.Circuit.AllBusNames()[k2]
                in_lines = self.circuit.lines.get_line_list(bus_name, 'in')  # upstream buses
                out_lines = self.circuit.lines.get_line_list(bus_name, 'out')  # downstream buses
                for cplx in range(0,2):
                    idxbs = dss.Circuit.AllBusNames().index(dss.Circuit.AllBusNames()[k2])      
                    if cplx == 0:               
                        load_val = load_kw[ph,idxbs]
                        cap_val = 0
                    else:  
                        load_val = load_kvar[ph,idxbs] 
                        cap_val = caparr[ph][idxbs]  
                    #gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))]) #some derivatives
                    hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
                                        [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])
                    available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                    if available_phases[ph] == 1:                 #quadratic terms
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2] = \
                                                                                                        -load_val * (beta_Z + (0.5 * beta_I* hessian_mag[0][0])) + \
                                                                                                        cap_val * (gamma_Z + (0.5 * gamma_I * hessian_mag[0][0]))# TE replace assignment w/ -load_val * beta_Z; #a**2               
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2 + 1] = \
                                                                                                        -load_val * (beta_Z + (0.5 * beta_I * hessian_mag[1][1])) + \
                                                                                                        cap_val * (gamma_Z + (0.5 * gamma_I * hessian_mag[1][1]))# TE replace assignment w/ -load_val * beta_Z; #b**2
                            #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1] = -load_val * beta_I * hessian_mag[0][1] / 2 #remove for TE
                            #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2] =  -load_val * beta_I * hessian_mag[1][0] / 2 #remove for TE

                    for i in range(len(in_lines)): # in lines            
                        line_idx = self.circuit.lines.get_idx(in_lines[i])
                        if available_phases[ph] == 1:
                            if cplx == 0: #real residual
                                #A_m and C_lm
                                H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx] = 1/2
                                H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2] = 1/2
                                #B_m and D_lm
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = 1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = 1/2
                            if cplx == 1: #imaginary residual
                                # #A_m, D_lm
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = -1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = -1/2
                                #B_m and C_lm
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx] = 1/2
                                H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = 1/2
                
                    for j in range(len(out_lines)): # out lines    
                        line_idx = self.circuit.lines.get_idx(out_lines[j])
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
        tf_no = len(dss.Transformers.AllNames()) - len(dss.RegControls.AllNames())
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
        vr_count = len(dss.RegControls.AllNames())
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
            for k2 in range(1, len(dss.Circuit.AllBusNames())):
                for cplx in range(0,2):
                    available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                    idxbs = dss.Circuit.AllBusNames().index(dss.Circuit.AllBusNames()[k2])
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
                    else:
                        if capacitance != 0 or der.imag != 0:
                            b_factor = capacitance - der.imag
                            b_factor = 0
                        else:                    
                            b_factor = caparr[ph][k2]                    
                    

                    if available_phases[ph] == 0: #if phase does not exist
                        b_temp = 0
                    else:
                        b_temp = (-load_val * beta_S) + (b_factor * gamma_S)
                    b[2*(nnode-1)*ph + 2*(k2-1) + cplx][0][0] = b_temp

            self.H, self.g, self.b = H, g, b

    def _init_KVL_matrices_vregs(self):
        """
        set H_reg, G_reg
        copied from 20180601/PYTHON/lib/compute_SBKVL_matrices.py
        written by @kathleenchang
        """
        # ---------- Voltage Regulator -----------
        tf_bus = self.circuit.transformers.get_bus_ph_matrix()
        vr_bus = self.circuit.voltage_regulators.get_bus_ph_matrix()
        tf_lines = self.circuit.transformers.get_num_lines_x_ph()
        tf_no = self.circuit.transformers.num_elements
        vr_lines = self.circuit.voltage_regulators.get_num_lines_x_ph()
        vr_count = self.circuit.voltage_regulators.num_elements
        nnode = self.circuit.buses.num_elements
        nline = self.circuit.lines.num_elements
        gain = self.circuit.voltage_regulators.get_gain_matrix()

        H_reg = np.zeros((2*vr_lines, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
        G_reg = np.zeros((2*vr_lines, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))

        #  voltage ratio: V_bus2 - gamma V_bus1 = 0
        line_in_idx = range(0, 2*vr_lines, 2)
        line_out_idx = range(1, 2*vr_lines, 2)

        vr_counter = 0
        for m in range(vr_count):
            bus1_idx = vr_bus[0, m]
            bus2_idx = vr_bus[1, m]
            for ph in range(0,3):
                if vr_bus[ph + 2,m] != 0:
                    # voltage gain: gamma*A_in = A_out
                    # gamma * B_in = B_out
                    # Negative not shown below because inserted into gain term
                    G_reg[2*vr_counter][2*nnode*ph + 2*bus1_idx] = gain[m] #A_in
                    G_reg[2*vr_counter][2*nnode*ph + 2*bus2_idx] = 1 #A_out
                    G_reg[2*vr_counter + 1][2*nnode*ph + 2*bus1_idx + 1] = gain[m]  #B_in
                    G_reg[2*vr_counter + 1][2*nnode*ph + 2*bus2_idx + 1] = 1 #B_out

                    #conservation of power: V_bus1 (I_bus1,out)* -  V_bus2 (I_bus2,in)* = 0
                    # A_1 * C_out + B_1 * D_out  - (A_2 * C_in + B_2 * D_in)
                    # j(B_1 * C_out - A_1 * D_out) - j(B_2 * C_in - A_2 * D_in)

                    #A_1 C_out
                    H_reg[2*vr_counter][2*nnode*ph + 2*bus1_idx][2*3*(nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter]] = 1
                    H_reg[2*vr_counter][2*3*(nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter]][2*nnode*ph + 2*bus1_idx] = 1
                    #B_1 D_out
                    H_reg[2*vr_counter][2*nnode*ph + 2*bus1_idx + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter] + 1] = 1
                    H_reg[2*vr_counter][2*3*(nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter]+ 1][2*nnode*ph + 2*bus1_idx + 1] = 1
                    #A_2 C_in
                    H_reg[2*vr_counter][2*nnode*ph + 2*bus2_idx][2*3*(nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter]] = -1
                    H_reg[2*vr_counter][2*3*(nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter]][2*nnode*ph + 2*bus2_idx] = -1
                    #B_2 D_in
                    H_reg[2*vr_counter][2*nnode*ph + 2*bus2_idx + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter] + 1] = -1
                    H_reg[2*vr_counter][2*3*(nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter] + 1][2*nnode*ph + 2*bus2_idx + 1] = -1

                    #B_1 * C_out
                    H_reg[2*vr_counter + 1][2*nnode*ph + 2*bus1_idx + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter]] = 1
                    H_reg[2*vr_counter + 1][2*3*(nnode+nline) + 2*tf_lines+ 2*line_out_idx[vr_counter]][2*nnode*ph + 2*bus1_idx + 1] = 1
                    # A_1 * D_out
                    H_reg[2*vr_counter + 1][2*nnode*ph + 2*bus1_idx][2*3*(nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter] + 1] = -1
                    H_reg[2*vr_counter + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_out_idx[vr_counter] + 1][2*nnode*ph + 2*bus1_idx] = -1
                    #B_2 * C_in
                    H_reg[2*vr_counter + 1][2*nnode*ph + 2*bus2_idx + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter]] = -1
                    H_reg[2*vr_counter + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter]][2*nnode*ph + 2*bus2_idx + 1] = -1
                    # A_2 * D_in
                    H_reg[2*vr_counter + 1][2*nnode*ph + 2*bus2_idx][2*3*(nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter] + 1] = 1
                    H_reg[2*vr_counter + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_in_idx[vr_counter] + 1][2*nnode*ph + 2*bus2_idx] = 1
                    vr_counter += 1
        self.H_reg, self.G_reg = H_reg, G_reg

    def compute_NR3FT(self):
        """ From 20180601/PYTHON/lib/compute_NR3FT.py, written by @kathleenchang """
        nnode = self.circuit.buses.num_elements
        nline = self.circuit.lines.num_elements
        vr_lines = self.circuit.voltage_regulators.get_num_lines_x_ph
        X = self.XNR
        g_SB, b_SB = self.g_SB, self.b_SB
        G_KVL, b_KVL = self.G_KVL, self.b_KVL
        H, g, b = self.H, self.g, self.b
        H_reg, G_reg = self.H_reg, self.G_reg

        FTSUBV = (g_SB @ X) - b_SB
        FTKVL = (G_KVL @ X) - b_KVL

        FTKCL = np.zeros((2*3*(nnode-1), 1))
        for i in range(2*3*(nnode-1)):
            r = (X.T @ (H[i, :, :] @ X)) \
            + (g[i, 0,:] @ X) \
            + b[i, 0,0]
            FTKCL[i, :] = r

        FTVR = np.zeros((2*vr_lines, 1)) #need to fix
        for i in range(2*vr_lines):
            r = X.T @ (H_reg[i, :, :] @ X)
            FTVR[i, :] = r

        FTVR2 = G_reg @ X

        FT = np.r_[FTSUBV, FTKVL, FTKCL, FTVR, FTVR2]

        return FT

    def compute_NR3JT(self):
        """ From 20180601/PYTHON/lib/compute_NR3JT.py, written by @kathleenchang """
        nnode = self.circuit.buses.num_elements
        nline = self.circuit.lines.num_elements
        vr_lines = self.circuit.voltage_regulators.get_num_lines_x_ph()
        tf_lines = self.circuit.transformers.get_num_lines_x_ph()

        X = self.XNR
        g_SB, b_SB = self.g_SB, self.b_SB
        G_KVL, b_KVL = self.G_KVL, self.b_KVL
        H, g, b = self.H, self.g, self.b
        H_reg, G_reg = self.H_reg, self.G_reg

        JSUBV = g_SB
        JKVL = G_KVL

        JKCL = np.zeros((2*3*(nnode-1), 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
        for i in range(2*3*(nnode-1)):
            r = (2 * (X.T @ H[i, :, :])) \
            + (g[i, 0, :])
            JKCL[i,:] = r

        JVR = np.zeros((2*vr_lines, 2*3*(nnode+nline)+ 2*tf_lines + 2*2*vr_lines))
        for i in range(2*vr_lines):
            r = (2 * (X.T @ H_reg[i, :, :]))
            JVR[i, :] = r

        JVR2 = G_reg

        JT = np.r_[JSUBV, JKVL, JKCL, JVR, JVR2]

        return JT

    def change_KCL_matrices(self, H, g, b, t, der, capacitance, wpu):

        nnode = self.circuit.buses.num_elements

        # load_kw, load_kvar = nominal_load_values(t)
        # caparr = nominal_cap_arr()

        load_kw = self.circuit.loads.get_ppu_matrix()
        load_kvar = self.circuit.loads.get_qpu_matrix()
        caparr = self.circuit.capacitors.get_cappu_matrix()

        # ----------Residuals for KCL at a bus (m) ----------

        #Zip Parameters
        beta_S = Circuit.aPQ_p
        beta_I = Circuit.aI_p
        beta_Z = Circuit.aZ_p

        # Capacitors
        gamma_S = Circuit.aPQ_p
        gamma_I = Circuit.aI_p
        gamma_Z = Circuit.aZ_p

        # Quadratic Terms

        # Time varying load
        if t != -1:
            for ph in range(0, 3):
                for k2 in range(1, nnode):  # skip slack bus
                    bus = self.circuit.buses.get_element(k2)
                    #dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k2]) #set the bus
                    available_phases = bus.phase_matrix  # phase array at specific bus
                    for cplx in range(0, 2):
                        if available_phases[ph] == 1:  # quadratic terms
                            H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2 *
                                                                                     (nnode)*ph + 2*k2] *= (1 + 0.1*np.sin(2*np.pi*0.01*t))
                            #-load_val * (beta_Z + (0.5 * beta_I* hessian_mag[0][0])) # TE replace assignment w/ -load_val * beta_Z; #a**2
                            H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2 *
                                                                                         (nnode)*ph + 2*k2 + 1] *= (1 + 0.1*np.sin(2*np.pi*0.01*t))
                            #-load_val * (beta_Z  + (0.5 * beta_I * hessian_mag[1][1])) # TE replace assignment w/ -load_val * beta_Z; #b**2

        # Constant Term
        if t != -1 or der != 0 or capacitance != 0:
            for ph in range(0, 3):
                for k2 in range(1, nnode):
                    bus = self.circuit.buses.get_element(k2)
                    available_phases = bus.phase_matrix
                    if available_phases[ph] == 1:
                        for cplx in range(0, 2):
                            if cplx == 0:
                                load_val = load_kw[ph][k2]
                                if der.real != 0:
                                    b_factor = der.real
                                else:
                                    b_factor = 0
                            else:
                                load_val = load_kvar[ph][k2]
                                if capacitance != 0 or der.imag != 0:
                                    b_factor = capacitance - der.imag
                                else:
                                    b_factor = caparr[ph][k2]
                            if t != -1:
                                print('Should not enter')
                                load_val *= (1 + 0.1*np.sin(2*np.pi*0.01*t))

                            b_temp = (-load_val * beta_S) + \
                                (b_factor * gamma_S)

                            b[2*(nnode-1)*ph + 2*(k2-1) +
                              cplx][0][0] = b_temp - wpu[ph][k2]

        return H, b

    def solve(self, time_delta, der=0, capacitance=0):
        """ From 20180601/PYTHON/lib/NR3.py, written by @Kathleen Chang"""
        dss, vvc_objects = self.dss, self.volt_var_controllers

        nnode = self.circuit.buses.num_elements
        nline = self.circuit.lines.num_elements
        vr_lines = self.circuit.voltage_regulators.get_num_lines_x_ph
        tf_lines = self.circuit.transformers.get_num_lines_x_ph

        XNR = self.XNR
        g_SB, b_SB = self.g_SB, self.b_SB
        G_KVL, b_KVL = self.G_KVL, self.b_KVL
        H, g, b = self.H, self.g, self.b
        H_reg, G_reg = self.H_reg, self.G_reg

        V0, I0 = self.__class__.V0, self.__class__.I0
        Vslack = self.__class__.VSLACK
        tol = self.__class__.tolerance
        maxiter = self.__class__.maxiter
        slacknode = self.__class__.SLACKIDX

        if maxiter == None:
            maxiter = 100

        FT = 1e99
        itercount = 0


        # adjust KCL based on capacitance, DER, and time-varying load
        # if der != 0 or capacitance != 0 or time_delta != -1:
        #     H, b = change_KCL_matrices(fn, H, g, b, time_delta, der, capacitance)

        # solve power-flow

        while np.amax(np.abs(FT)) >= 1e-9 and itercount < maxiter:
            print("Iteration number %f" % (itercount))
            FT = self.compute_NR3FT(XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode, nline, H_reg, G_reg, vr_lines)
            JT = self.compute_NR3JT(XNR, g_SB, G_KVL, H, g, nnode, nline, H_reg, G_reg, tf_lines, vr_lines)

            if JT.shape[0] >= JT.shape[1]:
                XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
            itercount += 1
        XNR_final = XNR

        XNR_compare = np.zeros(XNR_final.shape)
        ################ What a mess ....
        while np.linalg.norm(XNR_final - XNR_compare) > 1e-6:
            print('WPU')
            wpu = np.zeros((3, nnode))
            for vvo in vvc_objects:
                busName = vvo.get_busName()
                idxbus = dss.Circuit.AllBusNames().index(busName)
                phase = vvo.get_phase()
                qpu = vvo.get_Q(np.abs(XNR[2*nnode*phase + 2 * idxbus] + (1j*XNR[2*nnode*phase + 2 * idxbus+1])))
                print(qpu)
                wpu[phase, idxbus] = qpu

            H, b = self.change_KCL_matrices(H, g, b, time_delta, der, capacitance, wpu)

            while np.amax(np.abs(FT)) >= 1e-6 and itercount < maxiter:
                print("Iteration number %f" % (itercount))
                FT = self.compute_NR3FT(XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode, nline, H_reg, G_reg, vr_lines)
                JT = self.compute_NR3JT(XNR, g_SB, G_KVL, H, g, nnode, nline, H_reg, G_reg, tf_lines, vr_lines)

                if JT.shape[0] >= JT.shape[1]:
                    XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
                itercount += 1
            XNR_compare = XNR
        XNR_final = XNR_compare
        ##############

        # returns associated indices (values as list) of a bus's voltage regulators (keys)
        vr_idx_dict = self.circuit.voltage_regulators.voltage_regulator_index_dict
        vr_line_idx = range(0, vr_lines)

        # flag if need to rerun NR3
        flag = 0
        vr_line_counter = 0

        for k in vr_idx_dict.keys():
            for vridx in vr_idx_dict[k]: #{bus: [indices in dss.RegControls.AllNames(), ...]}

                dss.RegControls.Name(dss.RegControls.AllNames()[vridx])
                dss.Circuit.SetActiveBus(dss.CktElement.BusNames()[0].split(".")[0])
                winding = dss.RegControls.Winding()

                Vbase = dss.Bus.kVBase() *1000
                Sbase = 10**6
                Ibase = Sbase / Vbase
                band = dss.RegControls.ForwardBand()
                target_voltage = dss.RegControls.ForwardVreg()

                idxbs = dss.Circuit.AllBusNames().index((dss.CktElement.BusNames()[0].split('.')[0]))

                ph = dss.CktElement.BusNames()[0].split('.')[1:]
                ph_arr = [0, 0, 0]
                for i in ph:
                    ph_arr[int(i) - 1] = 1
                if len(ph) == 0:
                    ph_arr = [1, 1, 1]
                for ph in range(len(ph_arr)):
                    if ph_arr[ph] == 1: #loop over existing phases of voltage regulator

                        NR_voltage = np.abs((XNR[2*nnode*ph + 2*idxbs]  + (1j * XNR[2*nnode*ph + 2*idxbs + 1])) * Vbase / dss.RegControls.PTRatio())
                        print(NR_voltage)
                        if dss.RegControls.ForwardR() and dss.RegControls.ForwardX() and dss.RegControls.CTPrimary():
                            # if LDC exists

                            #vr_line_counter - counts the number of lines passed; two lines for every phase
                            #vridx - index of current voltage regulator in dss.RegControls.AllNames()
                            #tf_lines - number of transformers

                            line_idx =  2*vr_line_idx[vr_line_counter] + 2*(winding - 1)

                            I_reg = XNR[2*3*(nnode+nline) + 2*tf_lines + line_idx] + \
                                1j * XNR[2*3*(nnode+nline) + 2*tf_lines +  line_idx + 1]

                            V_drop = (dss.RegControls.ForwardR() + 1j*dss.RegControls.ForwardX()) / 0.2 * (I_reg * Ibase / dss.RegControls.CTPrimary())
                            #print(V_drop)
                            V_drop = (dss.RegControls.ForwardR() + 1j*dss.RegControls.ForwardX()) / 0.2 * (I_reg * Ibase / (dss.RegControls.CTPrimary()/0.2))
                            V_R = np.abs(NR_voltage - V_drop)

                            abs_diff = np.abs(V_R - target_voltage)
                            V_compare = V_R

                            print('V_R', V_R)
                            print('vdrop', V_drop)

                        else:
                            # if LDC term does not exist
                            print('**** LDC missing term ***** ')
                            abs_diff = np.abs(NR_voltage - target_voltage)
                            V_compare = NR_voltage

                        print('absolute difference: ', abs_diff, "\n")
                        vr_line_counter += 1

                        # compare NR3 voltage to forward Vreg voltage +- band
                        if abs_diff <= band: #converges
                            XNR_final = XNR
                            continue

                        elif abs_diff > band:
                            if V_compare > (target_voltage + band): #NR3 voltage above forward-Vreg
                                if dss.RegControls.TapNumber() <= -16 :
                                    print('Tap Number Out of Bounds' )
                                    XNR_final = XNR

                                else:
                                    print('Decrease Tap Number')
                                    dss.RegControls.TapNumber(dss.RegControls.TapNumber() - 1)
                                    print('New tap number ', dss.RegControls.TapNumber())
                                    flag = 1 #run NR3 again
                            else: #NR3 voltage below forward-Vreg
                                if dss.RegControls.TapNumber() >= 16:
                                    print('Tap Number Out of Bounds' )
                                    print('New tap number ', dss.RegControls.TapNumber())
                                    XNR_final = XNR

                                else:
                                    print('Increase tap number')
                                    dss.RegControls.TapNumber(dss.RegControls.TapNumber() + 1)
                                    flag = 1 #run NR3 again
            if flag == 1:
                print('\n Next iteration: ')
                XNR_final = self.solve(time_delta, der=0, capacitance=0)
                flag = 0

        self.XNR = XNR_final
