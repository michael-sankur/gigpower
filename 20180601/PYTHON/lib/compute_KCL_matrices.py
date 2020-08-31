import numpy as np
import opendssdirect as dss
import re
import sys
def compute_KCL_matrices(fn, t, der, capacitance):

    dss.run_command('Redirect ' + fn)
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())
    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])
    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0

    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase

    def load_order_f():
        load_order = {}
        for n in range(len(dss.Loads.AllNames())):
            dss.Loads.Name(dss.Loads.AllNames()[n])
            pattern =  r"(\w+)\."
            load_bus = re.findall(pattern, dss.CktElement.BusNames()[0])
            if load_bus[0] not in load_order:
                load_order[load_bus[0]] = 1
            elif load_bus[0] in load_order:
                load_order[load_bus[0]] += 1
        return load_order
    load_order_list = load_order_f()

    def load_values():
        load_ph_arr = np.zeros((nnode, max(load_order_list.values()), 3))
        load_kw_arr_ph = np.zeros((3, nnode))
        load_kvar_arr_ph = np.zeros((3, nnode))
        if t == -1:
            var = 1
        else:
            var = (1 + 0.1*np.sin(2*np.pi*0.01*t))

        for load in range(len(dss.Loads.AllNames())):
            dss.Loads.Name(dss.Loads.AllNames()[load])
            pattern =  r"(\w+)\."
            load_bus = re.findall(pattern, dss.CktElement.BusNames()[0])
            load_ph_arr_temp = [0, 0, 0]
            for i in range(1, 4):
                pattern = r"\.%s" % (str(i))
                load_ph = re.findall(pattern, dss.CktElement.BusNames()[0])
                if load_ph:
                    load_ph_arr_temp[i - 1] = 1
            for j in range(max(load_order_list.values())):
                idxbs = dss.Circuit.AllBusNames().index(load_bus[0])
                if np.all(load_ph_arr[idxbs, j,:] == [0, 0, 0]):
                    load_ph_arr[idxbs, j, :] = load_ph_arr_temp
                    for i in range(len(load_ph_arr_temp)):
                        if load_ph_arr_temp[i] == 1:
                            load_kw_arr_ph[i,idxbs] += dss.Loads.kW() *1e3*var / Sbase / sum(load_ph_arr_temp)
                            load_kvar_arr_ph[i,idxbs] += dss.Loads.kvar()*1e3*var / Sbase  / sum(load_ph_arr_temp)
                    break
        return load_kw_arr_ph, load_kvar_arr_ph
    load_kw_arr_ph, load_kvar_arr_ph = load_values()

    def cap_dict():
        cap_dict = {}
        cap_ph_dict = {}
        caparr = np.zeros((3, nnode))
        for n in range(len(dss.Capacitors.AllNames())):
            dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
            load_phases = [0, 0, 0]
            for i in range(1, 4): #if the phase is present, what other phases are
                pattern = r"\.%s" % (str(i))
                m2 = re.findall(pattern, dss.CktElement.BusNames()[0])
                if m2:
                    load_phases[i - 1] = 1
            cap_ph_dict[dss.CktElement.BusNames()[0]] = load_phases
            cap_dict[dss.CktElement.BusNames()[0]] = dss.Capacitors.kvar() * 1e3 / Sbase / sum(load_phases)
        for n in cap_ph_dict.keys():
            pattern =  r"(\w+)\."
            cap_bus = re.findall(pattern, n)
            idxbs = dss.Circuit.AllBusNames().index(cap_bus[0])
            for ph in range(3):
                if cap_ph_dict[n][ph] == 1:
                    caparr[ph, idxbs] += cap_dict[n]
        return cap_dict, cap_ph_dict,caparr

    cap_dict_kvar, cap_dict_kv, caparr = cap_dict()

    def bus_phases():
        dictionary = {}
        for k2 in range(len(dss.Circuit.AllNodeNames())):
            for i in range(1, 4):
                pattern = r"\.%s" % (str(i))
                m = re.findall(pattern, dss.Circuit.AllNodeNames()[k2])
                a, b = dss.Circuit.AllNodeNames()[k2].split('.')
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

    def get_line_idx(line): #returns the index of a line as stored in dss.Lines.AllNames()
        k = -1
        for n in range(len(dss.Lines.AllNames())):
            if dss.Lines.AllNames()[n] == line:
                k = n
        return k

    def linelist(busname): #returns two lists of in and out lines at a bus
        in_lines = np.array([])
        out_lines = np.array([])
        for k in range(len(dss.Lines.AllNames())):
            dss.Lines.Name(dss.Lines.AllNames()[k])
            if busname in dss.Lines.Bus1():
                out_lines = np.append(out_lines, dss.Lines.AllNames()[k])
            elif busname in dss.Lines.Bus2():
                in_lines = np.append(in_lines, dss.Lines.AllNames()[k])
        return in_lines,out_lines

    # ----------Residuals for KCL at a bus (m) ----------
    bp = bus_phases()

    beta_S = 1
    beta_I = 0.0
    beta_Z = 0.0

    H = np.zeros((2*3*(nnode-1), 2 * 3 * (nnode + nline), 2 * 3* (nnode + nline)))
    g = np.zeros((2*3*(nnode-1), 1, 2*3*(nnode+nline)))
    b = np.zeros((2*3*(nnode-1), 1, 1))

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
            in_lines, out_lines = linelist(dss.Circuit.AllBusNames()[k2]) #get in/out lines of bus
            for cplx in range(0,2):
                idxbs = dss.Circuit.AllBusNames().index(dss.Circuit.AllBusNames()[0])
                if cplx == 0:
                    load_val = load_kw_arr_ph[ph][idxbs]
                else:
                    load_val = load_kvar_arr_ph[ph][idxbs]
                gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))]) #some derivatives
                hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
                                    [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])
                available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                if available_phases[ph] == 1:                 #quadratic terms
                    H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2] = -load_val * (beta_Z)# + (0.5 * beta_I* hessian_mag[0][0])) # TE replace assignment w/ -load_val * beta_Z; #a**2
                    H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2 + 1] = -load_val * (beta_Z)#  + (0.5 * beta_I * hessian_mag[1][1])) # TE replace assignment w/ -load_val * beta_Z; #b**2
                    #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1] = -load_val * beta_I * hessian_mag[0][1] / 2 #remove for TE
                    #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2] =  -load_val * beta_I * hessian_mag[1][0] / 2 #remove for TE

                for i in range(len(in_lines)): #fill in H for the inlines
                    line_idx = get_line_idx(in_lines[i])
                    if available_phases[ph] == 1:
                        if cplx == 0: #real residual
                            #A_m and C_lm
                            H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx] = 1/2
                            H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2] = 1/2
                            #B_m and D_lm
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = 1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = 1/2
                        if cplx == 1: #complex residual
                            #A_m, D_lm
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = -1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = -1/2
                            #B_m and C_lm
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx] = 1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = 1/2

                for j in range(len(out_lines)): #fill in H for the outlines
                    line_idx = get_line_idx(out_lines[j])
                    if available_phases[ph] == 1:
                        if cplx == 0:
                            #A_m and C_mn
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx] = -1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2] = -1/2
                            #B_m and D_mn
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = -1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = -1/2
                        if cplx == 1:
                            #A_m and D_mn
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1]= 1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = 1/2
                            #C_m and B_mn
                            H[2*ph*(nnode-1) + (k2-1)*2+cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx] = -1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = -1/2

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
                #load_val = d_factor(dss.Circuit.AllBusNames()[k2], cplx, ph)
                idxbs = dss.Circuit.AllBusNames().index(dss.Circuit.AllBusNames()[0])
                if cplx == 0:
                    load_val = load_kw_arr_ph[ph][idxbs]
                else:
                    load_val = load_kvar_arr_ph[ph][idxbs]
                #linear terms
                g_temp = np.zeros(2*3*(nnode+nline)) #preallocate g
                available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                if available_phases[ph] == 0: #if phase does not exist
                    g_temp[2*(ph)*nnode + 2*k2 + cplx] = 1
                else:
                    h = 5
                    #g_temp[2*ph*nnode+ 2 * k2] = -load_val * beta_I * ((-(A0 * hessian_mag[0][0] + B0 * hessian_mag[0][1])) \
                                          # +  gradient_mag[0]) #remove for TE
                    #g_temp[2*ph*nnode+ 2 * k2 + 1] = -load_val * beta_I * ((-(A0 * hessian_mag[0][1] + B0 * hessian_mag[1][1])) \
                                        #       +  gradient_mag[1]) #remove for TE
                g[2*(nnode-1)*ph + 2*(k2-1) + cplx, 0,:] = g_temp

                #constant terms
                b_factor = 0 #DER term
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
                        b_factor = caparr[ph][k2]
                else:
                    b_factor = 0


                if available_phases[ph] == 0: #if phase does not exist at bus, set b = 0
                    b_temp = 0
                else:
                    b_temp = (-load_val * beta_S) + b_factor #TE version
                    # b_temp = -load_val * (beta_S \
                    # + (beta_I) * (((hessian_mag[0][1] * A0 * B0) + ((1/2)*hessian_mag[0][0] * ((A0)**2)) + ((1/2)*hessian_mag[1][1] * (B0**2))) \
                    # -  (A0 * gradient_mag[0] + B0* gradient_mag[1]) \
                    # +  (A0**2 + B0**2) ** (1/2))) \
                    # + b_factor #calculate out the constant term in the residual

                b[2*(nnode-1)*ph + 2*(k2-1) + cplx][0][0] = b_temp #store the in the b matrix

    return H, g, b
