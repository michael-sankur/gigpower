import numpy as np
import opendssdirect as dss
from lib.compute_NR3FT_vectorized import compute_NR3FT_vectorized as ft1
from lib.compute_NR3JT_vectorized import compute_NR3JT_vectorized as jt1
from lib.compute_vecmat import compute_vecmat
from lib.relevant_openDSS_parameters import relevant_openDSS_parameters
import re

#a31faff

class LinDist3Flow:
    def __init__(self, dss_path, slacknode, Vslack, V0=None, I0=None, tol=1e-9, maxiter=100):

        dss.run_command('Redirect ' + '"' + dss_path + '"')
        dss.Solution.Solve()
        nline = len(dss.Lines.AllNames())
        nnode = len(dss.Circuit.AllBusNames())

        self.all_line_names = dss.Lines.AllNames()
        self.nline = len(self.all_line_names)

        self.all_bus_names = dss.Circuit.AllBusNames()
        self.nnode = len(self.all_bus_names)

        self.all_node_names = dss.Circuit.AllNodeNames()
        self.nn = len(self.all_node_names)

        XNR = np.zeros((2*3*(nnode + nline),1))


        # intialize node voltage portion of XNR
        if V0 == None or len(V0) == 0:
            for ph in range(0,3):
                for k1 in range(0,nnode):

                    XNR[2*ph*nnode + 2*k1] = Vslack[ph].real
                    XNR[2*ph*nnode + 2*k1+1] = Vslack[ph].imag

        # If initial V is given (usually from CVX)
        elif len(V0) != 0:
            for ph in range(0,3):
                for k1 in range(0,nnode):
                    XNR[2*ph*nnode + 2*k1] = V0[ph,k1].real
                    XNR[2*ph*nnode + 2*k1+1] = V0[ph,k1].imag


        # intialize line current portion of XNR
        if I0 == None or len(I0) == 0:
            for k1 in range(0,nnode):
                XNR[(2*3*nnode):] = 0.0*np.ones((6*nline,1))

        elif len(I0) != 0:
            for ph in range(0,3):
                for k1 in range(0,nline):
                    XNR[(2*3*nnode) + 2*ph*nline + 2*k1] = I0[ph,k1].real
                    XNR[(2*3*nnode) + 2*ph*nline + 2*k1+1] = I0[ph,k1].imag

        # getting additional info from openDSS
        self.max_iter = maxiter

        # additional info from openDSS
        XNR1, g_SB, b_SB, G_KVL, b_KVL, H, g, b = compute_vecmat(XNR, dss_path, Vslack)
        TXnum, RXnum, PH, spu, APQ, AZ, AI, cappu, wpu, vvcpu = relevant_openDSS_parameters(dss_path)

        self.XNR = XNR
        self.g_SB = g_SB # constant
        self.b_SB = b_SB # constant
        self.G_KVL = G_KVL # constant
        self.b_KVL = b_KVL # constant
        self.H = H
        self.g = g
        self.b = b

        self.TXnum = TXnum # constant
        self.RXnum = RXnum # constant
        self.PH = PH # constant
        self.spu = spu # constant
        self.APQ = APQ # constant
        self.AZ = AZ # constant
        self.AI = AI # constant
        self.cappu = cappu # constant
        self.wpu = wpu # constant
        self.vvcpu = vvcpu # constant

        self._is_init = True

        self.bus_load = np.zeros((self.nnode, 3, 2))  #bus_name, phase, kwkvar

        # in_lines, out_lines
        self.in_lines = []
        self.out_lines = []

        for k2 in range(1, self.nnode):
            in_lines, out_lines = self._linelist(self.all_bus_names[k2])
            self.in_lines.append(in_lines)
            self.out_lines.append(out_lines)
        self.bus_phases = self._bus_phases()

    def solve(self):

        if self._is_init:
            self._is_init = False
        else:
            self.H, self.g, self.b = self._calculate_Hgb()

        # TODO: update XNR, H, g, b after solve()
        iter_count = 0 #count # of iterations
        XNR_deep_vec = np.zeros((self.max_iter+1, 2*3*(self.nnode+self.nline), 1))
        FT_deep_vec = np.zeros((self.max_iter+1,2*3*(self.nnode+self.nline), 1))
        JT_deep_vec = np.zeros((self.max_iter+1,2*3*(self.nnode+self.nline), 2*3*(self.nnode+self.nline)))
        FTSUBV_vec = np.zeros((self.max_iter+1,6, 1))
        FTKVL_vec = np.zeros((self.max_iter+1,36, 1))
        FTKCL_vec = np.zeros((self.max_iter+1,36, 1))
        JSUBV_vec = np.zeros((self.max_iter+1,6, 78))
        JKVL_vec = np.zeros((self.max_iter+1,36, 78))
        JKCL_vec = np.zeros((self.max_iter+1,36, 78))

        XNR_deep_vec[iter_count, :, 0] = self.XNR[:, 0]
        FT1 = 1e99
        while np.amax(np.abs(FT1)) >= 1e-9 and iter_count < self.max_iter:
            FT1 = ft1(self.XNR, self.g_SB, self.b_SB, self.G_KVL, self.b_KVL, self.H, self.g, self.b, self.nnode) #vectorized
            FT_deep_vec[iter_count,:, :] = FT1
            FTSUBV_vec[iter_count, :, 0] = np.reshape(FT1[0:6], (6))
            FTKVL_vec[iter_count,:, 0] = np.reshape(FT1[6:42], (36))
            FTKCL_vec[iter_count, :, 0] = np.reshape(FT1[42:], (36))

            JT1 = jt1(self.XNR, self.g_SB, self.G_KVL, self.H, self.g, self.nnode, self.nline)
            JT_deep_vec[iter_count,:, :] = JT1
            JSUBV_vec[iter_count,:, :] = np.reshape(JT1[0:6], (6, 78))
            JKVL_vec[iter_count,:, :] = np.reshape(JT1[6:42], (36, 78))
            JKCL_vec[iter_count,:, :] = np.reshape(JT1[42:], (36, 78))

            iter_count += 1 #count up at this point (indexing should match)

            # if JT.shape[0] >= JT.shape[1]: #non-vectorized
            #     XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
            if JT1.shape[0] >= JT1.shape[1]: #vectorized
                self.XNR = self.XNR - np.linalg.inv(JT1.T@JT1)@JT1.T@FT1

            #dump xnr into the mega-XNR matrix
            XNR_deep_vec[iter_count, :, 0] = self.XNR[:, 0]

        ###############################################################################
        VNR = np.zeros((3,self.nnode), dtype='complex')
        for ph in range(0,3):
            for k1 in range(0,self.nnode):
                VNR[ph,k1] = self.XNR[2*ph*self.nnode + 2*k1] + 1j*self.XNR[2*ph*self.nnode + 2*k1+1]
                if np.abs(VNR[ph,k1].real) <= 1e-12:
                    VNR[ph,k1] = 0 + VNR[ph,k1].imag
                if np.abs(VNR[ph,k1].imag) <= 1e-12:
                    VNR[ph,k1] = VNR[ph,k1].real + 0
        VNR[self.PH == 0] = 0
        XNR = self.XNR[2*3*self.nnode:]
        INR = np.zeros((3,self.nline), dtype='complex')
        for ph in range(0,3):
            for k1 in range(0,self.nline):
                INR[ph,k1] = XNR[2*ph*self.nline + 2*k1] + 1j*XNR[2*ph*self.nline + 2*k1+1]
                if np.abs(INR[ph,k1].real) <= 1e-12:
                    INR[ph,k1] = 0 + INR[ph,k1].imag
                if np.abs(INR[ph,k1].imag) <= 1e-12:
                    INR[ph,k1] = INR[ph,k1].real + 0

        STXNR = np.zeros((3,self.nnode), dtype='complex')
        SRXNR = np.zeros((3,self.nnode), dtype='complex')
        for ph in range(0,3):
            for k1 in range(0,self.nline):
                STXNR[ph,k1] = VNR[ph,self.TXnum[k1]]*np.conj(INR[ph,k1])
                if np.abs(STXNR[ph,k1].real) <= 1e-12:
                    STXNR[ph,k1] = 0 + STXNR[ph,k1].imag
                if np.abs(STXNR[ph,k1].imag) <= 1e-12:
                    STXNR[ph,k1] = STXNR[ph,k1].real + 0
                SRXNR[ph,k1] = VNR[ph,self.RXnum[k1]]*np.conj(INR[ph,k1]) #needs to be updated
                if np.abs(SRXNR[ph,k1].real) <= 1e-12:
                    SRXNR[ph,k1] = 0 + SRXNR[ph,k1].imag
                if np.abs(SRXNR[ph,k1].imag) <= 1e-12:
                    SRXNR[ph,k1] = SRXNR[ph,k1].real + 0

        sNR = np.zeros((3,self.nnode), dtype='complex')
        iNR = np.zeros((3,self.nnode), dtype='complex')
        # Totdal node loads
        sNR = self.spu*(self.APQ + self.AI*np.abs(VNR) + self.AZ*np.abs(VNR)**2) - 1j*self.cappu.real + self.wpu + 1j*self.vvcpu.real
        sNR[self.PH == 0] = 0
        for ph in range(0,3):
            for k1 in range(0,self.nnode):
                if np.abs(sNR[ph,k1].real) <= 1e-12:
                    sNR[ph,k1] = 0 + sNR[ph,k1].imag
                if np.abs(sNR[ph,k1].imag) <= 1e-12:
                    sNR[ph,k1] = sNR[ph,k1].real + 0
        sNR[self.PH == 0] = 0
        # Total node current

        iNR[self.PH != 0] = np.conj(sNR[self.PH != 0]/VNR[self.PH != 0]) #also needs to be updated...
        iNR[self.PH == 0] = 0
        for ph in range(0,3):
            for k1 in range(0,self.nnode):
                if np.abs(iNR[ph,k1].real) <= 1e-12:
                    iNR[ph,k1] = 0 + iNR[ph,k1].imag
                if np.abs(iNR[ph,k1].imag) <= 1e-12:
                    iNR[ph,k1] = iNR[ph,k1].real + 0
        iNR[self.PH == 0] = 0

        #reset bus_load, ready for the next timestep
        self.bus_load = np.zeros((self.nnode, 3, 2))  #bus_name, phase, kwkvar
        # VNR, INR, iNR, sNR, STXNR, SRXNR

    def set_bus_kw(self, bus, kw, ph):
        bus_idx = self.all_bus_names.index(bus)
        self.bus_load[bus_idx, ph, 0] = kw

    def set_bus_kvar(self, bus, kvar, ph):
        bus_idx = self.all_bus_names.index(bus)
        self.bus_load[bus_idx, ph, 1] = kvar

    def _calculate_Hgb(self):
        beta_S = 0.9
        beta_I = 0
        beta_Z = 0.1
        H = np.zeros((2*3*(self.nnode-1), 2 * 3 * (self.nnode + self.nline), 2 * 3* (self.nnode + self.nline)))
        g = np.zeros((2*3*(self.nnode-1), 1, 2*3*(self.nnode+self.nline)))
        b = np.zeros((2*3*(self.nnode-1), 1, 1))

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
            for k2 in range(1, self.nnode): #skip slack bus
                in_lines = self.in_lines[k2-1]
                out_lines = self.out_lines[k2-1]

                for cplx in range(0,2):
                    load_val = self.bus_load[k2, ph, cplx]
                    gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))]) #some derivatives
                    hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
                                        [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])

                    available_phases = self.bus_phases[self.all_bus_names[k2]] #phase array at specific bus
                    if available_phases[ph] == 1:                 #quadratic terms
                        H[2*ph*(self.nnode-1) + (k2-1)*2 + cplx][2*(self.nnode)*ph + 2*k2][2*(self.nnode)*ph + 2*k2] = -load_val * (beta_Z) #+ (0.5 * beta_I* hessian_mag[0][0])) # TE replace right side of equality with -load_val * beta_Z #a**2
                        H[2*ph*(self.nnode-1) + (k2-1)*2 + cplx][2*(self.nnode)*ph + 2*k2 + 1][2*(self.nnode)*ph + 2*k2 + 1] = -load_val * (beta_Z)# + (0.5 * beta_I * hessian_mag[1][1])) # TE -load_val * beta_Z #b**2
                        #H[2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1][2*ph*(nnode-1) + (k2-1)*2] = -load_val * beta_I * hessian_mag[0][1] #cross quad. terms in taylor exp,TE  remove
                        #H[2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2][2*ph*(nnode-1) + (k2-1)*2] =  -load_val * beta_I * hessian_mag[0][1] #TE remove

                    for i in range(len(in_lines)): #fill in H for the inlines
                        line_idx = self.all_line_names.index(in_lines[i])
                        if available_phases[ph] == 1:
                            if cplx == 0: #real residual
                                #A_m and C_lm
                                H[2*ph*(self.nnode-1) + (k2-1)*2 + cplx][2*(self.nnode)*ph + 2*k2][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx] = 1/2
                                H[2*ph*(self.nnode-1) + (k2-1)*2 + cplx][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx][2*(self.nnode)*ph + 2*k2] = 1/2
                                #B_m and D_lm
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*(self.nnode)*ph + 2*k2 + 1][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx + 1] = 1/2
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx + 1][2*(self.nnode)*ph + 2*k2 + 1] = 1/2
                            if cplx == 1: #complex residual
                                #A_m, D_lm
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*(self.nnode)*ph + 2*k2][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx + 1] = -1/2
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx + 1][2*(self.nnode)*ph + 2*k2] = -1/2
                                #B_m and C_lm
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*(self.nnode)*ph + 2*k2 + 1][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx] = 1/2
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx][2*(self.nnode)*ph + 2*k2 + 1] = 1/2

                    for j in range(len(out_lines)): #fill in H for the outlines
                        line_idx = self.all_line_names.index(out_lines[j])

                        if available_phases[ph] == 1:
                            if cplx == 0:
                                #A_m and C_mn
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*(self.nnode)*ph + 2*k2][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx] = -1/2
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx][2*(self.nnode)*ph + 2*k2] = -1/2
                                #B_m and D_mn
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*(self.nnode)*ph + 2*k2 + 1][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx + 1] = -1/2
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx + 1][2*(self.nnode)*ph + 2*k2 + 1] = -1/2
                            if cplx == 1:
                                #A_m and D_mn
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*(self.nnode)*ph + 2*k2][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx + 1]= 1/2
                                H[2*ph*(self.nnode-1) + (k2-1)*2+ cplx][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx + 1][2*(self.nnode)*ph + 2*k2] = 1/2
                                #C_m and B_mn
                                H[2*ph*(self.nnode-1) + (k2-1)*2+cplx][2*(self.nnode)*ph + 2*k2 + 1][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx] = -1/2
                                H[2*ph*(self.nnode-1) + (k2-1)*2+cplx][2*3*(self.nnode) + 2*ph*self.nline + 2*line_idx][2*(self.nnode)*ph + 2*k2 + 1] = -1/2

        #Linear Term
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
            for k2 in range(1, len(dss.Circuit.AllBusNames())):
                for cplx in range(0,2):
                    load_val = self.bus_load[k2, ph, cplx]
                    #linear terms
                    g_temp = np.zeros(len(self.XNR)) #preallocate g
                    available_phases = self.bus_phases[self.all_bus_names[k2]] #phase array at specific bus

                    if available_phases[ph] == 0: #if phase does not exist
                        g_temp[2*(ph)*self.nnode + 2*k2 + cplx] = 1
                        # g_temp[2*ph*nnode+ 2 * k2] = -load_val* beta_I * ((1/2 * (-2 * A0 * hessian_mag[0][0] - 2 * B0 * hessian_mag[0][1])) #remove lines for TE\
                        #                       +  gradient_mag[0]) # TE remove
                        # g_temp[2*ph*nnode+ 2 * k2 + 1] = -load_val * beta_I * ((1/2 * (-2* A0 *hessian_mag[0][1] - 2 * B0 * hessian_mag[1][1])) #remove lines \
                        #                           +  gradient_mag[1]) #TE remove
                    g[2*(self.nnode-1)*ph + 2*(k2-1) + cplx, 0,:] = g_temp #o.w.

                    #constant terms
                    b_factor = 0
                    if cplx == 0:
                        b_factor = 0 #DER term
                    elif cplx == 1:
                        #b_factor = (dss.Capacitors.kvar()*1000/Sbase) #DER term
                        b_factor = 0

                    if available_phases[ph] == 0: #if phase does not exist at bus, set b = 0
                        b_temp = 0
                    else:
                        b_temp = (-load_val * beta_S) + b_factor

                    #-load_val * (beta_S\ /
                        # + (beta_I) * (hessian_mag[0][1] * A0 * B0 + (1/2)*hessian_mag[0][0] * ((A0)**2) + (1/2)*hessian_mag[1][1] * (B0**2)) \
                        # - beta_I * (A0 * gradient_mag[0] +B0* gradient_mag[1]) \
                        # + beta_I * (A0**2 + B0**2) ** (1/2)) \
                        # + b_factor #calculate out the constant term in the residual
                        ##-load_val * beta_S + b_factor #TE version
                    b[2*(self.nnode-1)*ph + 2*(k2-1) + cplx][0][0] = b_temp #store the in the b matrix

        return H, g, b

    def _bus_phases(self): #goes through all the buses and saves their phases to a list stored in a dictionary
        #1 if phase exists, 0 o.w.
        #list goes [a, b, c]
        #key is the bus name (without the phase part)
        dictionary = {}
        for k2 in range(self.nn):
            for i in range(1, 4):
                pattern = r"\.%s" % (str(i))

                m = re.findall(pattern, self.all_node_names[k2])
                a, b = self.all_node_names[k2].split('.')
                if m and a in dictionary:
                    temp = dictionary[a]
                    temp[i - 1] = 1
                    dictionary[a] = temp
                elif m and a not in dictionary:
                    dictionary[a] = [0, 0, 0]
                    temp = dictionary[a]
                    temp[i - 1] = 1
                    dictionary[a] = temp
        return dictionary

    def _linelist(self, busname): #returns two lists of in and out lines at a bus
        dss.Circuit.SetActiveBus(busname)
        in_lines = np.array([])
        out_lines = np.array([])
        for k in range(self.nline):
            dss.Lines.Name(self.all_line_names[k])
            if busname in dss.Lines.Bus1():
                out_lines = np.append(out_lines, self.all_line_names[k])
            elif busname in dss.Lines.Bus2():
                in_lines = np.append(in_lines, self.all_line_names[k])
        return in_lines, out_lines