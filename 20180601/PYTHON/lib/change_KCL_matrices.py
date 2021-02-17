import numpy as np
import opendssdirect as dss
import re
import sys
def change_KCL_matrices(fn, H, g, b, t, der, capacitance):

    #dss.run_command('Redirect ' + fn)
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())
    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])
    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0

    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase

    # 3 x nnode array of loads
    def load_values():
        load_kw_arr_ph = np.zeros((3, nnode))
        load_kvar_arr_ph = np.zeros((3, nnode))
        if t == -1:
            var = 1
        else:
            var = (1 + 0.1*np.sin(2*np.pi*0.01*t))
        for load in range(len(dss.Loads.AllNames())):
            dss.Loads.Name(dss.Loads.AllNames()[load])
            load_data = dss.CktElement.BusNames()[0].split('.')
            idxbs = dss.Circuit.AllBusNames().index(load_data[0])
            for ph in range(1, len(load_data)):
                load_kw_arr_ph[int(load_data[ph]) - 1, idxbs] += dss.Loads.kW() * 1e3 *var / Sbase / (len(load_data) - 1)
                load_kvar_arr_ph[int(load_data[ph]) - 1, idxbs] += dss.Loads.kvar() * 1e3 * var / Sbase / (len(load_data) - 1)
        return load_kw_arr_ph, load_kvar_arr_ph
    load_kw_arr_ph, load_kvar_arr_ph = load_values()

    # 3 x nnode array of capacitance
    def cap_arr():
        caparr = np.zeros((3, nnode))
        for n in range(len(dss.Capacitors.AllNames())):
            dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
            cap_data = dss.CktElement.BusNames()[0].split('.')

            idxbs = dss.Circuit.AllBusNames().index(cap_data[0])
            for ph in range(1, len(cap_data)):
                caparr[int(cap_data[ph]) - 1, idxbs] += dss.Capacitors.kvar() * 1e3 / Sbase / (len(cap_data) - 1)
        return caparr

    caparr = cap_arr()

    # {bus : [1 x 3 phase existence]}
    def bus_phases():
        dictionary = {}
        for k2 in range(len(dss.Circuit.AllNodeNames())):
            a, b = dss.Circuit.AllNodeNames()[k2].split('.')
            if a in dictionary:
                temp = dictionary[a]
                temp[int(b) - 1] = 1
                dictionary[a] = temp
            elif a not in dictionary:
                dictionary[a] = [0, 0, 0]
                temp = dictionary[a]
                temp[int(b) - 1] = 1
                dictionary[a] = temp
        return dictionary

    # ----------Residuals for KCL at a bus (m) ----------
    bp = bus_phases()

    #Zip Parameters
    beta_S = 1.0
    beta_I = 0.0
    beta_Z = 0.0

    # Quadratic Terms
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
                available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                for cplx in range(0,2):
                    if available_phases[ph] == 1: #quadratic terms
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2] *= (1 + 0.1*np.sin(2*np.pi*0.01*t))#-load_val * (beta_Z + (0.5 * beta_I* hessian_mag[0][0])) # TE replace assignment w/ -load_val * beta_Z; #a**2
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2 + 1]*= (1 + 0.1*np.sin(2*np.pi*0.01*t))#-load_val * (beta_Z  + (0.5 * beta_I * hessian_mag[1][1])) # TE replace assignment w/ -load_val * beta_Z; #b**2
                        #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1] = -load_val * beta_I * hessian_mag[0][1] / 2 #remove for TE
                        #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2] =  -load_val * beta_I * hessian_mag[1][0] / 2 #remove for TE

    # Constant Term
    if t != -1 or der != 0 or capacitance != 0:
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
                available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                idxbs = dss.Circuit.AllBusNames().index(dss.Circuit.AllBusNames()[k2])
                b_factor = 0
                for cplx in range(0,2):
                    if cplx == 0:
                        load_val = load_kw_arr_ph[ph][idxbs]
                        if der.real != 0:
                            b_factor = der.real
                        else:
                            b_factor = 0
                    else:
                        load_val = load_kvar_arr_ph[ph][idxbs]
                        if capacitance != 0 or der.imag != 0:
                            b_factor = capacitance - der.imag
                        else:
                            b_factor = caparr[ph][k2]

                    if available_phases[ph] == 0: #if phase does not exist at bus, set b = 0
                        b_temp = 0
                    else:
                        b_temp = (-load_val * beta_S) + b_factor #TE version

                    b[2*(nnode-1)*ph + 2*(k2-1) + cplx][0][0] = b_temp #store the in the b matrix

    return H, b
