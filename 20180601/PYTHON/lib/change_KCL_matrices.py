import numpy as np
import opendssdirect as dss
import re
import sys
def change_KCL_matrices(fn, H, g, b, t, der, capacitance):

    dss.run_command('Redirect ' + fn)
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())
    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])
    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0

    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase

    def d_factor(busname, cplx, ph): #pass in bus phase, re/im, phase
        #hsould work for multiphase
        for n in range(len(dss.Loads.AllNames())): #go through the loads
            dss.Loads.Name(dss.Loads.AllNames()[n]) #set the load
            if busname in dss.CktElement.BusNames()[0]: #check is the busname in the busname of the load
                pattern =  r"\.%s" % (str(ph + 1)) #if it is, is the phase present?
                m = re.findall(pattern, dss.CktElement.BusNames()[0])
                if m:
                    load_phases = [0, 0, 0]
                    for i in range(1, 4): #if the phase is present, what other phases are
                        pattern = r"\.%s" % (str(i))
                        m = re.findall(pattern, dss.CktElement.BusNames()[0])
                        if m:
                            load_phases[i - 1] = 1
                    if t == -1:
                        var = 1
                    else:
                        var = (1 + 0.1*np.sin(2*np.pi*0.01*t))#adjust this later

                    if cplx == 0: #figure it out for multiphase
                        if sum(load_phases) == 1:
                            return dss.Loads.kW()*1*var*1e3/Sbase
                        else:
                            return (dss.Loads.kW())*1*var*1e3/Sbase/sum(load_phases)
                    else:
                        if sum(load_phases) == 1:
                            return dss.Loads.kvar()*1*var*1e3/Sbase
                        else:
                            return (dss.Loads.kvar())*1*var*1e3/Sbase/sum(load_phases)

        return 0

    def cap_dict():
        cap_dict_kvar = {}
        cap_dict_kv = {}
        for n in range(len(dss.Capacitors.AllNames())):
            dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
            #print(dss.CktElement.BusNames()[0])
            pattern =  r"(\w+)\."  #if it is, is the phase present?
            m = re.findall(pattern, dss.CktElement.BusNames()[0])
            cap_dict_kvar[m[0]] = dss.Capacitors.kvar()
            cap_dict_kv[m[0]] = dss.Capacitors.kV()
            #print(dss.Capacitors.kvar()*1000/Sbase)
        return cap_dict_kvar, cap_dict_kv
        #print(cap_dict)
    cap_dict_kvar, cap_dict_kv = cap_dict()

    def bus_phases(): #goes through all the buses and saves their phases to a list stored in a dictionary
    #1 if phase exists, 0 o.w.
    #list goes [a, b, c]
    #key is the bus name (without the phase part)
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

    # ----------Residuals for KCL at a bus (m) ----------
    bp = bus_phases()

    beta_S = 1
    beta_I = 0.0
    beta_Z = 0.0

    if t != -1:
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
                for cplx in range(0,2):
                    load_val = d_factor(dss.Circuit.AllBusNames()[k2], cplx, ph)
                    gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))]) #some derivatives
                    hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
                                        [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])
                    available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                    if available_phases[ph] == 1:                 #quadratic terms
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2] = -load_val * (beta_Z + (0.5 * beta_I* hessian_mag[0][0])) # TE replace assignment w/ -load_val * beta_Z; #a**2
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2 + 1] = -load_val * (beta_Z  + (0.5 * beta_I * hessian_mag[1][1])) # TE replace assignment w/ -load_val * beta_Z; #b**2
                        #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1] = -load_val * beta_I * hessian_mag[0][1] / 2 #remove for TE
                        #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2] =  -load_val * beta_I * hessian_mag[1][0] / 2 #remove for TE

    if t!= -1 or der != 0 or capacitance != 0:
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
                    if t != -1:
                        load_val = d_factor(dss.Circuit.AllBusNames()[k2], cplx, ph)
                        #linear terms
                        g_temp = np.zeros(2*3*(nnode+nline)) #preallocate g
                        available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                        if available_phases[ph] == 0: #if phase does not exist
                            g_temp[2*(ph)*nnode + 2*k2 + cplx] = 1
                        else:
                            h = 6
                            # g_temp[2*ph*nnode+ 2 * k2] = -load_val * beta_I * ((-(A0 * hessian_mag[0][0] + B0 * hessian_mag[0][1])) \
                            #                        +  gradient_mag[0]) #remove for TE
                            # g_temp[2*ph*nnode+ 2 * k2 + 1] = -load_val * beta_I * ((-(A0 * hessian_mag[0][1] + B0 * hessian_mag[1][1])) \
                            #                            +  gradient_mag[1]) #remove for TE
                        g[2*(nnode-1)*ph + 2*(k2-1) + cplx, 0,:] = g_temp

                    #constant terms

                    b_factor = 0 #DER term
                    if cplx == 0:
                        #add line about u der.real
                        if der.real != 0:
                            b_factor = der.real
                        else:
                            b_factor = 0
                    elif cplx == 1:
                        #add line about capacitance & reactive portion (v) der.imag
                        if capacitance != 0 or der.imag != 0:
                            b_factor = capacitance - der.imag
                        else:
                            if dss.Circuit.AllBusNames()[k2] in cap_dict_kvar.keys():
                                if dss.Circuit.AllBusNames()[k2] == '675':
                                    print('679 classic song')
                                    b_factor = cap_dict_kvar[dss.Circuit.AllBusNames()[k2]] * 1000 / Sbase /3
                                else:
                                    b_factor = cap_dict_kvar[dss.Circuit.AllBusNames()[k2]] * 1000 / Sbase

                            else:
                                #b_factor = (dss.Capacitors.kvar()*1000/Sbase) #DER term
                                b_factor = 0
                        print("\n")
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
