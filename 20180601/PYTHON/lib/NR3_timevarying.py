import numpy as np
import opendssdirect as dss
import time
import re

from lib.compute_NR3FT_vectorized import compute_NR3FT_vectorized as ft
from lib.compute_NR3JT_vectorized import compute_NR3JT_vectorized as jt
from lib.change_KCL_matrices import change_KCL_matrices
from lib.helper import transformer_regulator_parameters, voltage_regulator_index_dict
from lib.map_output import map_output
from lib.basematrices import basematrices

def NR3_timevarying(fn, slacknode, Vslack, V0, I0, tol, maxiter, der, capacitance, time_delta):
    XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg = basematrices(fn, slacknode, Vslack, V0, I0)
   #need to move time, cap, DER into this function
    nline = len(dss.Lines.AllNames())  
    nnode = len(dss.Circuit.AllBusNames())

    if tol == None:
        tol = 1e-9

    if maxiter == None:
        maxiter = 100

    FT = 1e99
    itercount = 0

    _, _, tf_lines, vr_lines, _, _, _ = transformer_regulator_parameters()
    
    # adjust KCL based on capacitance, DER, and time-varying load
    # if der != 0 or capacitance != 0 or time_delta != -1:  
    #     H, b = change_KCL_matrices(fn, H, g, b, time_delta, der, capacitance)

    # solve power-flow
    while np.amax(np.abs(FT)) >= 1e-9 and itercount < maxiter:
        print("Iteration number %f" % (itercount))
        FT = ft(XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode, nline, H_reg, G_reg, vr_lines)
        JT = jt(XNR, g_SB, G_KVL, H, g, nnode, nline, H_reg, G_reg, tf_lines, vr_lines)
    
        if JT.shape[0] >= JT.shape[1]:
            XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
        itercount+=1
 
    
    # vr_idx_dict = voltage_regulator_index_dict()
    # flag = 0
    # for k in vr_idx_dict.keys():
    #     for p in vr_idx_dict[k]:
            
    #         dss.RegControls.Name(dss.RegControls.AllNames()[p])
    #         #dss.Transformers.Name(dss.RegControls.AllNames()[p])
    #         #won't necessarily have the same name
    #         print(dss.RegControls.Name())
    #         tw = dss.RegControls.TapWinding()
    #         #dss.Circuit.SetActiveBus(dss.CktElement.BusNames()[tw-1].split(".")[0])
    #         dss.Circuit.SetActiveBus(dss.CktElement.BusNames()[0].split(".")[0])
        
    #         Vbase = dss.Bus.kVBase() * 1000
    #         band = dss.RegControls.ForwardBand()
    #         target_voltage = dss.RegControls.ForwardVreg()

    #         #idxbs = dss.Circuit.AllBusNames().index((dss.CktElement.BusNames()[tw - 1].split('.')[0]))
    #         idxbs = dss.Circuit.AllBusNames().index((dss.CktElement.BusNames()[0].split('.')[0]))
        
    #         #ph = dss.CktElement.BusNames()[tw - 1].split('.')[1:] 
    #         ph = dss.CktElement.BusNames()[0].split('.')[1:] 
    #         ph_arr = [0, 0, 0]
    #         for i in ph:
    #             ph_arr[int(i) - 1] = 1
    #         if len(ph) == 0:
    #             ph_arr = [1, 1, 1]

    #         XNR_final = 0
    #         for ph in range(len(ph_arr)):
    #             #print(ph)
    #             if ph_arr[ph] == 1:
    #                 NR_voltage = np.abs(XNR[2*nnode*ph + 2*idxbs] + 1j*XNR[2*nnode*ph + 2*idxbs + 1]) * Vbase / dss.RegControls.PTRatio()
    #                 print(NR_voltage)

    #                 diff = np.abs(NR_voltage - target_voltage)
    #                 print(diff)

    #                 # print('NR voltage: ', NR_voltage)
    #                 # print('difference: ', diff)
    #                 # print(ph, idxbs, dss.Bus.Name(), band, target_voltage)

    #                 #compare voltage to target voltage
    #                 #assess the difference
    #                 if diff <= band: #converges
    #                     print('Converges')
    #                     XNR_final = XNR
    #                     continue
                            
    #                 elif diff > band:
    #                     if NR_voltage > (target_voltage + band): 
    #                         if dss.RegControls.TapNumber() <= -16 :
    #                             print('Tap Number Out of Bounds' )
    #                             XNR_final = XNR
                        
    #                         else:
    #                             print('decrease Tap Number')
    #                             dss.RegControls.TapNumber(dss.RegControls.TapNumber() - 1)
    #                             print('new tapnumber ', dss.RegControls.TapNumber())
    #                             flag = 1
    #                     else:
    #                         if dss.RegControls.TapNumber() >= 16:
    #                             print('Tap Number Out of Bounds' )
    #                             print('new tapnumber ', dss.RegControls.TapNumber())
    #                             XNR_final = XNR
                        
    #                         else:
    #                             print('increase tap number')
    #                             dss.RegControls.TapNumber(dss.RegControls.TapNumber() + 1)
    #                             flag = 1
    #     if flag == 1:  
    #         print('different iteration ')                
    #         XNR_final = NR3_timevarying(fn, slacknode, Vslack, V0, I0, tol, maxiter, der, capacitance, time_delta)
    #         flag = 0

    return XNR
   