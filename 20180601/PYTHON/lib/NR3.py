import numpy as np
import opendssdirect as dss
import time
import re

from lib.compute_NR3FT import compute_NR3FT 
from lib.compute_NR3JT import compute_NR3JT 
from lib.change_KCL_matrices import change_KCL_matrices
from lib.helper import transformer_regulator_parameters, voltage_regulator_index_dict
from lib.map_output import map_output
from lib.basematrices import basematrices

def NR3(fn, slacknode, Vslack, V0, I0, tol, maxiter, der, capacitance, time_delta, vvc_objects):
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
        FT = compute_NR3FT(XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode, nline, H_reg, G_reg, vr_lines)
        JT = compute_NR3JT(XNR, g_SB, G_KVL, H, g, nnode, nline, H_reg, G_reg, tf_lines, vr_lines)
    
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

        H, b = change_KCL_matrices(H, g, b, time_delta, der, capacitance, wpu)

        while np.amax(np.abs(FT)) >= 1e-6 and itercount < maxiter:
            print("Iteration number %f" % (itercount))
            FT = compute_NR3FT(XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode, nline, H_reg, G_reg, vr_lines)
            JT = compute_NR3JT(XNR, g_SB, G_KVL, H, g, nnode, nline, H_reg, G_reg, tf_lines, vr_lines)
        
            if JT.shape[0] >= JT.shape[1]:
                XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
            itercount += 1
        XNR_compare = XNR
    XNR_final = XNR_compare
    ##############

    # returns associated indices (values as list) of a bus's voltage regulators (keys)
    vr_idx_dict = voltage_regulator_index_dict() 
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
            XNR_final = NR3(fn, slacknode, Vslack, V0, I0, tol, maxiter, der, capacitance, time_delta, vvc_objects)
            flag = 0

    return XNR_final
   