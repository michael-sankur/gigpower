import numpy as np
from lib.compute_vecmatY import compute_vecmat
from lib.compute_KCL_matricesY import compute_KCL_matrices
import opendssdirect as dss
import re
import time
def basematrices(fn, slacknode, Vslack, V0, I0):

    dss.run_command('Redirect ' + fn)

    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())
    #t0 = time.time()
    XNR = np.zeros((2*3*(nnode),1))

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

    g_SB, b_SB, G_KVL, b_KVL, R_matrix, X_matrix, G_matrix, Hb_matrix = compute_vecmat(fn, Vslack)

    H, g, b = compute_KCL_matrices(fn, -1, 0, 0,  R_matrix, X_matrix, G_matrix, Hb_matrix)

    return XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b
