import numpy as np
from lib.compute_vecmat import compute_vecmat
from lib.compute_KCL_matrices import compute_KCL_matrices
from lib.helper import transformer_regulator_parameters
import opendssdirect as dss
import re
import time
def basematrices(fn, slacknode, Vslack, V0, I0):
    tf_bus, vr_bus, tf_lines, vr_lines, tf_count, vr_no, gain = transformer_regulator_parameters()

    nline = len(dss.Lines.AllNames())  
    nnode = len(dss.Circuit.AllBusNames())
    XNR = np.zeros((2*3*(nnode + nline) + 2*tf_lines + 2*2*vr_lines,1))


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
        XNR[(2*3*nnode):] = 0.0*np.ones((6*nline + 2*tf_lines + 2*2*vr_lines ,1))

    # If initial I is given
    elif len(I0) != 0:
        for ph in range(0,3):
            for k1 in range(0,len(dss.Lines.AllNames())):
                XNR[(2*3*nnode) + 2*ph*nline + 2*k1] = I0[ph,k1].real
                XNR[(2*3*nnode) + 2*ph*nline + 2*k1+1] = I0[ph,k1].imag
        XNR[(2*3*nnode + 2*3*nline):] = np.zeros((len(XNR) - 2*3*nnode - 2*3*nline), 1)

    # generate static matrices
    XNR, g_SB, b_SB, G_KVL, b_KVL, H_reg, G_reg = compute_vecmat(XNR, fn, Vslack, tf_bus, vr_bus, tf_lines, vr_lines, tf_count, vr_no, gain)
    # generate non-static matrices
    H, g, b = compute_KCL_matrices(fn, -1, 0, 0, tf_bus, vr_bus, tf_lines, vr_lines)
    return XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg
