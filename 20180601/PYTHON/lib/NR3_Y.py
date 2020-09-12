import numpy as np
from lib.compute_NR3FT_vectorized0 import compute_NR3FT_vectorized as ft1
from lib.compute_NR3JT_vectorized0 import compute_NR3JT_vectorized as jt1
from lib.relevant_openDSS_parameters import relevant_openDSS_parameters
import opendssdirect as dss
import time
import re
def NR3_timevarying(fn, XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, tol, maxiter, der, capacitance, time_delta):
    dss.run_command('Redirect ' + fn)
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())

    if tol == None:
        tol = 1e-9

    if maxiter == None:
        maxiter = 100

    FT = 1e99
    itercount = 0
    # solve power-flow
    while np.amax(np.abs(FT)) >= 1e-9 and itercount < maxiter:
        print("Iteration number %f" % (itercount))
        FT = ft1(XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode)
        JT = jt1(XNR, g_SB, G_KVL, H, g, nnode, nline)

        if JT.shape[0] >= JT.shape[1]:
            XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
        itercount+=1
    # # retrieve relevant parameters for formatting output
    TXnum, RXnum, PH, spu, APQ, AZ, AI, cappu, wpu, vvcpu = \
        relevant_openDSS_parameters(fn, time_delta)

    #remap XNR to VNR, INR, STXNR, SRXNR, iNR, sNR
    #VNR = XNR(1:2:2*3*nnode-1).' + 1j*XNR(2:2:2*3*nnode).';
    VNR = np.zeros((3,nnode), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nnode):
            VNR[ph,k1] = XNR[2*ph*nnode + 2*k1] + 1j*XNR[2*ph*nnode + 2*k1+1]
            if np.abs(VNR[ph,k1].real) <= 1e-12:
                VNR[ph,k1] = 0 + VNR[ph,k1].imag
            if np.abs(VNR[ph,k1].imag) <= 1e-12:
                VNR[ph,k1] = VNR[ph,k1].real + 0
    VNR[PH == 0] = 0
    XNR = XNR[2*3*nnode:]
    print('VNR')
    print(VNR)

    return VNR
