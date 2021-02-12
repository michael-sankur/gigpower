import numpy as np
import opendssdirect as dss
import time
import re

from lib.compute_NR3FT_vectorized import compute_NR3FT_vectorized as ft
from lib.compute_NR3JT_vectorized import compute_NR3JT_vectorized as jt
from lib.change_KCL_matrices import change_KCL_matrices
from lib.helper import transformer_regulator_parameters, simple_reg_control
from lib.map_output import map_output

def NR3_timevarying(fn, XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, tol, maxiter, der, capacitance, time_delta, H_reg, G_reg):
   # dss.run_command('Redirect ' + fn)
   
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
 
    XNR = simple_reg_control(XNR)

    VNR, INR, STXNR, SRXNR, iNR, sNR = map_output(nline, nnode, XNR, fn, time_delta)

    return VNR, INR, STXNR, SRXNR, iNR, sNR, itercount
