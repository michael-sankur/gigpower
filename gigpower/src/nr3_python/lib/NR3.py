import numpy as np
import opendssdirect as dss
import time
import re
import pandas as pd

from . compute_NR3FT import compute_NR3FT 
from . compute_NR3JT import compute_NR3JT 
from . change_KCL_matrices import change_KCL_matrices
from . helper import transformer_regulator_parameters, wpu_final_arr
from . map_output import map_output
from . basematrices import basematrices
from . regulator_control import regulator_control

def NR3(slacknode, Vslack, V0, I0, tol, maxiter, der, capacitance, time_delta, vvc_objects):
    XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg = basematrices(slacknode, Vslack, V0, I0)
   #need to move time, cap, DER into this function
    nline = len(dss.Lines.AllNames())  
    nnode = len(dss.Circuit.AllBusNames())
    Sbase = 1e6

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
        print("Iteration number for Original NR3 %f" % (itercount))
        FT = compute_NR3FT(XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode, nline, H_reg, G_reg, vr_lines)
        JT = compute_NR3JT(XNR, g_SB, G_KVL, H, g, nnode, nline, H_reg, G_reg, tf_lines, vr_lines)
    
        if JT.shape[0] >= JT.shape[1]:
            XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
        itercount += 1
    XNR_final = XNR

    # print non-VVC results
    VNR, INR, STXNR, SRXNR, iNR, sNR = map_output(nline, nnode, XNR, 0, [])


    XNR_current = np.zeros(XNR_final.shape)
    XNR_previous = XNR_final
 
    ############ WPU / Volt VAR Control Iterations
    # itercount = 0
    # while np.linalg.norm(XNR_previous - XNR_current) > 1e-7 and itercount < 20: 
    #     print("Iteration number for VVC %f" % (itercount), np.linalg.norm(XNR_previous - XNR_current) )
               
    #     if itercount != 0: #if you have already iterated
    #         XNR_previous = XNR_current

    #     wpu = np.zeros((3, nnode))

    #     for vvo in vvc_objects: 
    #         # iterate through VVC objects 
    #         # and collect data to find reactive power for VVC
    #         busName = vvo.get_busName()
    #         idxbus = dss.Circuit.AllBusNames().index(busName)
    #         phase = vvo.get_phase()
    #         qpu = vvo.get_Q(np.abs(XNR_previous[2*nnode*phase + 2 * idxbus] + \
    #             (1j*XNR_previous[2*nnode*phase + 2 * idxbus+1])))

    #         wpu[phase, idxbus] = qpu*1e3/Sbase

    #     # adjust KCL matrices based on VVC
    #     H, b = change_KCL_matrices(H, g, b, time_delta, der, capacitance, wpu)

    #     # comparison variables
    #     XNR_in_loop = XNR_previous
    #     FT = 1e99

    #     while np.amax(np.abs(FT)) >= 1e-6:
            
    #         FT = compute_NR3FT(XNR_in_loop, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode, nline, H_reg, G_reg, vr_lines)
    #         JT = compute_NR3JT(XNR_in_loop, g_SB, G_KVL, H, g, nnode, nline, H_reg, G_reg, tf_lines, vr_lines)
        
    #         if JT.shape[0] >= JT.shape[1]:
    #             XNR_in_loop = XNR_in_loop - np.linalg.inv(JT.T@JT)@JT.T@FT
          
    #     XNR_current = XNR_in_loop
    #     itercount += 1
    
    # # WPU End
    # XNR_final = XNR_current
  
    ##############
    #Regulator Control
    # flag, XNR_final = regulator_control(XNR, vr_lines, tf_lines, nnode, nline)
    
    # if flag == 1:  
    #     print('\n Regulator Control next iteration: ')                
    #     XNR_final = NR3(slacknode, Vslack, V0, I0, tol, maxiter, der, capacitance, time_delta, vvc_objects)
    #     flag = 0
    
    ##############
    # print out VVC results
    # VNR2, INR, STXNR, SRXNR, iNR, sNR = map_output(nline, nnode, XNR_final, 0, vvc_objects)
    # print('VNR before VVC', np.abs(VNR))
    # print('VNR after VVC', np.abs(VNR2))
    # print('Absolute  VNR - VNR2', np.abs(VNR)-np.abs(VNR2))
    
    # # dataframe print out
    # VDSS_df = pd.DataFrame(data=np.abs(VNR),columns=[dss.Circuit.AllBusNames()],index=['A', 'B', 'C']).T
    # DSS_df = pd.DataFrame(data=np.abs(VNR2),columns=[dss.Circuit.AllBusNames()],index=['A', 'B', 'C']).T
    # print(VDSS_df)
    # print(DSS_df)
    # print(VDSS_df - DSS_df)

    return XNR_final
   