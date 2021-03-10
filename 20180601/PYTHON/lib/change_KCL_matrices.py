import numpy as np
import opendssdirect as dss
import re
import sys
from lib.helper import bus_phases, nominal_load_values, nominal_cap_arr
from lib.zip_parameters import *

def change_KCL_matrices(H, g, b, t, der, capacitance, wpu):

    nnode = len(dss.Circuit.AllBusNames())

    load_kw, load_kvar = nominal_load_values(t) 
    caparr = nominal_cap_arr()

    # ----------Residuals for KCL at a bus (m) ----------
    bp = bus_phases()

    #Zip Parameters
    beta_S = get_S()
    beta_I = get_I()
    beta_Z = get_Z()

    gamma_S = 0.0
    gamma_I = 0.0
    gamma_Z = 1

    # Quadratic Terms

    # Time varying load
    if t != -1:
        for ph in range(0,3):
            for k2 in range(1, len(dss.Circuit.AllBusNames())): #skip slack bus
                dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k2]) #set the bus
                available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                for cplx in range(0,2):
                    if available_phases[ph] == 1: #quadratic terms
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2] *= (1 + 0.1*np.sin(2*np.pi*0.01*t))
                        #-load_val * (beta_Z + (0.5 * beta_I* hessian_mag[0][0])) # TE replace assignment w/ -load_val * beta_Z; #a**2
                        H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2 + 1]*= (1 + 0.1*np.sin(2*np.pi*0.01*t))
                        #-load_val * (beta_Z  + (0.5 * beta_I * hessian_mag[1][1])) # TE replace assignment w/ -load_val * beta_Z; #b**2

    # wpu 
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
                idxbs = dss.Circuit.AllBusNames().index(dss.Circuit.AllBusNames()[k2])      
                if cplx == 0:               
                    load_val = load_kw[ph,idxbs]
                    cap_val = 0
                else:  
                    load_val = load_kvar[ph,idxbs]  
                    cap_val = caparr[ph][idxbs]  - wpu[ph][idxbs] #WPU HERE !!
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
                        

    # Constant Term
    if t != -1 or der != 0 or capacitance != 0:
        for ph in range(0,3):
            for k2 in range(1, len(dss.Circuit.AllBusNames())):
                available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                idxbs = dss.Circuit.AllBusNames().index(dss.Circuit.AllBusNames()[k2])
                for cplx in range(0,2):
                    if cplx == 0:
                        load_val = load_kw[ph][idxbs]
                        if der.real != 0:
                            b_factor = der.real
                        else:
                            b_factor = 0
                    else:
                        load_val = load_kvar[ph][idxbs]
                        if capacitance != 0 or der.imag != 0:
                            b_factor = capacitance - der.imag
                        else:
                            b_factor = caparr[ph][k2] - wpu[ph][k2]

                    if available_phases[ph] == 0: #if phase does not exist at bus, set b = 0
                        b_temp = 0
                    else:
                        b_temp = (-load_val * beta_S) + (b_factor * gamma_S) #TE version

                    b[2*(nnode-1)*ph + 2*(k2-1) + cplx][0][0] = b_temp #store the in the b matrix

    return H, b
