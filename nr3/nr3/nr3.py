import numpy as np
import opendssdirect as dss
from nr3.lib.compute_NR3FT_vectorized import compute_NR3FT_vectorized as ft1
from nr3.lib.compute_NR3JT_vectorized import compute_NR3JT_vectorized as jt1
from nr3.lib.compute_vecmat import compute_vecmat
from nr3.lib.relevant_openDSS_parameters import relevant_openDSS_parameters
import re

#a31faff
SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])

class NR3:
    def __init__(self, dss_path, slacknode=SLACKIDX, Vslack=VSLACK, V0=None, I0=None, tol=1e-9, maxiter=100):

        dss.run_command('Redirect ' + '"' + dss_path + '"')
        nline = len(dss.Lines.AllNames())
        nnode = len(dss.Circuit.AllBusNames())

        self.all_line_names = dss.Lines.AllNames()
        self.nline = len(self.all_line_names)

        self.all_bus_names = dss.Circuit.AllBusNames()
        self.nnode = len(self.all_bus_names)

        self.all_node_names = dss.Circuit.AllNodeNames()
        self.nn = len(self.all_node_names)

        self.all_load_names = dss.Loads.AllNames()
        self.nload = len(self.all_load_names)

        self.Vbase = dss.Bus.kVBase() * 1000
        self.Sbase = 1000000.0

        self.Ibase = self.Sbase/self.Vbase
        self.Zbase = self.Vbase/self.Ibase

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

        self.max_iter = maxiter

        # additional info from openDSS
        self.bus_load, self.bus_load_divide, self.load_ph_arr = self._init_load_values()  #bus_name, phase, kwkvar
        XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b = compute_vecmat(XNR, dss_path, Vslack, self.bus_load)
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
        self.spu = self.bus_load[:, :, 0] + 1j * self.bus_load[:, :, 1] # constant
        self.APQ = APQ # constant
        self.AZ = AZ # constant
        self.AI = AI # constant
        self.cappu = cappu # constant
        self.wpu = wpu # constant
        self.vvcpu = vvcpu # constant

        # in_lines, out_lines
        self.in_lines = []
        self.out_lines = []

        for k2 in range(1, self.nnode):
            in_lines, out_lines = self._linelist(self.all_bus_names[k2])
            self.in_lines.append(in_lines)
            self.out_lines.append(out_lines)
        self.bus_phases = self._bus_phases()

    def solve(self):

        self.H, self.g, self.b = self._calculate_Hgb()
        self.spu = self.bus_load[:, :, 0] + 1j * self.bus_load[:, :, 1] # constant
        # TODO: update XNR, H, g, b after solve()
        iter_count = 0 #count # of iterations

        FT1 = 1e99
        while np.amax(np.abs(FT1)) >= 1e-9 and iter_count < self.max_iter:
            FT1 = ft1(self.XNR, self.g_SB, self.b_SB, self.G_KVL, self.b_KVL, self.H, self.g, self.b, self.nnode) #vectorized
            JT1 = jt1(self.XNR, self.g_SB, self.G_KVL, self.H, self.g, self.nnode, self.nline)
            # if JT.shape[0] >= JT.shape[1]: #non-vectorized
            #     XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
            if JT1.shape[0] > JT1.shape[1]: #vectorized
                self.XNR = self.XNR - np.linalg.inv(JT1.T@JT1)@JT1.T@FT1

            iter_count += 1 #count up at this point (indexing should match)
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

        STXNR = np.zeros((3,self.nline), dtype='complex')
        SRXNR = np.zeros((3,self.nline  ), dtype='complex')
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
        # Total node current

        iNR[self.PH != 0] = np.conj(sNR[self.PH != 0]/VNR[self.PH != 0]) #also needs to be updated...
        iNR[self.PH == 0] = 0
        for ph in range(0,3):
            for k1 in range(0,self.nnode):
                if np.abs(iNR[ph,k1].real) <= 1e-12:
                    iNR[ph,k1] = 0 + iNR[ph,k1].imag
                if np.abs(iNR[ph,k1].imag) <= 1e-12:
                    iNR[ph,k1] = iNR[ph,k1].real + 0

        #reset bus_load, ready for the next timestep
        self.bus_load = np.zeros((3, self.nnode, 2))  #bus_name, phase, kwkvar

        return VNR, INR, iNR, sNR, STXNR, SRXNR

    def set_bus_kw(self, bus, kw, ph):
        bus_idx = self.all_bus_names.index(bus)
        self.bus_load[ph, bus_idx, 0] = kw/self.bus_load_divide[ph, bus_idx, 0]

    def set_bus_kvar(self, bus, kvar, ph):
        bus_idx = self.all_bus_names.index(bus)
        self.bus_load[ph, bus_idx, 1] = kvar/self.bus_load_divide[ph, bus_idx, 1]

    def set_load_kw(self, load, kw):
        load_idx = self.all_load_names.index(load)
        # TODO: refactor this
        try:
            bus = load.split('_')[1]
        except:
            bus = re.findall("\d+", load)[0]
        bus_idx = self.all_bus_names.index(bus)
        load_ph = self.load_ph_arr[load_idx]
        for i, ph in enumerate(load_ph):
            if ph == 1:
                self.bus_load[i, bus_idx, 0] += kw *1e3 / self.Sbase / sum(load_ph)

    def set_load_kvar(self, load, kvar):
        load_idx = self.all_load_names.index(load)
        # TODO: refactor this
        try:
            bus = load.split('_')[1]
        except:
            bus = re.findall("\d+", load)[0]
        bus_idx = self.all_bus_names.index(bus)
        load_ph = self.load_ph_arr[load_idx]
        for i, ph in enumerate(load_ph):
            if ph == 1:
                self.bus_load[i, bus_idx, 1] += kvar *1e3 / self.Sbase / sum(load_ph)

    def _calculate_Hgb(self):
        beta_S = 1.0
        beta_I = 0.0
        beta_Z = 0.0
        H = self.H
        g = self.g
        b = self.b

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
                for cplx in range(0,2):
                    idxbs = self.all_bus_names.index(self.all_bus_names[k2])
                    load_val = self.bus_load[ph, idxbs, cplx]
                    gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))]) #some derivatives
                    hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
                                        [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])

                    available_phases = self.bus_phases[self.all_bus_names[k2]] #phase array at specific bus
                    if available_phases[ph] == 1:                 #quadratic terms
                        H[2*ph*(self.nnode-1) + (k2-1)*2 + cplx][2*(self.nnode)*ph + 2*k2][2*(self.nnode)*ph + 2*k2] = -load_val * (beta_Z + (0.5 * beta_I* hessian_mag[0][0])) # TE replace assignment w/ -load_val * beta_Z; #a**2
                        H[2*ph*(self.nnode-1) + (k2-1)*2 + cplx][2*(self.nnode)*ph + 2*k2 + 1][2*(self.nnode)*ph + 2*k2 + 1] = -load_val * (beta_Z  + (0.5 * beta_I * hessian_mag[1][1])) # TE replace assignment w/ -load_val * beta_Z; #b**2

                        #H[2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1][2*ph*(nnode-1) + (k2-1)*2] = -load_val * beta_I * hessian_mag[0][1] #cross quad. terms in taylor exp,TE  remove
                        #H[2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2][2*ph*(nnode-1) + (k2-1)*2] =  -load_val * beta_I * hessian_mag[0][1] #TE remove
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
            for k2 in range(1, self.nnode):
                for cplx in range(0,2):
                    idxbs = self.all_bus_names.index(self.all_bus_names[k2])
                    load_val = self.bus_load[ph, idxbs, cplx]
                    #linear terms
                    g_temp = np.zeros(2*3*(self.nnode+self.nline)) #preallocate g
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


    def _init_load_order_f(self):
        load_order = {}
        for n in range(self.nload):
            dss.Loads.Name(self.all_load_names[n])
            pattern =  r"(\w+)\."
            load_bus = re.findall(pattern, dss.CktElement.BusNames()[0])
            if load_bus[0] not in load_order:
                load_order[load_bus[0]] = 1
            elif load_bus[0] in load_order:
                load_order[load_bus[0]] += 1
        return load_order


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