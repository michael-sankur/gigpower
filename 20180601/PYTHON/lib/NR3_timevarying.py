import numpy as np
#vectorized
from lib.compute_NR3FT_vectorized import compute_NR3FT_vectorized as ft1
from lib.compute_NR3JT_vectorized import compute_NR3JT_vectorized as jt1
from lib.compute_vecmat import compute_vecmat
from lib.change_KCL_matrices import change_KCL_matrices
from lib.relevant_openDSS_parameters import relevant_openDSS_parameters
import opendssdirect as dss
import re

def NR3_timevarying(fn, XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, tol, maxiter, der, capacitance, time):
    dss.run_command('Redirect ' + fn)
    #dss.Solution.Solve()
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())

    if tol == None:
        tol = 1e-9

    if maxiter == None:
        maxiter = 100

    FT1 = 1e99
    itercount = 0


    if der != 0 or capacitance != 0 or time != 1:
        H, g, b = change_KCL_matrices(fn, H, g, b, time, der, capacitance)
    while np.amax(np.abs(FT1)) >= 1e-9 and itercount < maxiter:
        print("Iteration number %f" % (itercount))
        # H, g, b = compute_KCL_matrices(fn, times[itercount], 0, 0)
        FT1 = ft1(XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode) #vectorized
        JT1 = jt1(XNR, g_SB, G_KVL, H, g, nnode, nline)

        if JT1.shape[0] >= JT1.shape[1]: #vectorized
            XNR = XNR - np.linalg.inv(JT1.T@JT1)@JT1.T@FT1

        itercount+=1 #diff from the other iter_count, limits # of iterations

    TXnum, RXnum, PH, spu, APQ, AZ, AI, cappu, wpu, vvcpu = \
        relevant_openDSS_parameters(fn)

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
    # INR = XNR(2*3*nnode+1:2:2*3*nnode+2*3*nline-1) + 1j*XNR(2*3*nnode+2:2:2*3*nnode+2*3*nline)
    INR = np.zeros((3,nline), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nline):
            INR[ph,k1] = XNR[2*ph*nline + 2*k1] + 1j*XNR[2*ph*nline + 2*k1+1]
            if np.abs(INR[ph,k1].real) <= 1e-12:
                INR[ph,k1] = 0 + INR[ph,k1].imag
            if np.abs(INR[ph,k1].imag) <= 1e-12:
                INR[ph,k1] = INR[ph,k1].real + 0
    # print('inr')
    print("INR:")
    print(INR)

    # STXNR_n^phi = V_m^phi (I_mn^phi)^*
    # SRXNR_n^phi = V_n^phi (I_mn^phi)^*
    STXNR = np.zeros((3,nline), dtype='complex')
    SRXNR = np.zeros((3,nline), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nline):
            STXNR[ph,k1] = VNR[ph,TXnum[k1]]*np.conj(INR[ph,k1])
            if np.abs(STXNR[ph,k1].real) <= 1e-12:
                STXNR[ph,k1] = 0 + STXNR[ph,k1].imag
            if np.abs(STXNR[ph,k1].imag) <= 1e-12:
                STXNR[ph,k1] = STXNR[ph,k1].real + 0
            SRXNR[ph,k1] = VNR[ph,RXnum[k1]]*np.conj(INR[ph,k1]) #needs to be updated
            if np.abs(SRXNR[ph,k1].real) <= 1e-12:
                SRXNR[ph,k1] = 0 + SRXNR[ph,k1].imag
            if np.abs(SRXNR[ph,k1].imag) <= 1e-12:
                SRXNR[ph,k1] = SRXNR[ph,k1].real + 0

    # print('stxnr and srxnr')
    print("STXNR:")
    print(STXNR)
    print("SRXNR:")
    print(SRXNR)
    # print("\n")


    sNR = np.zeros((3,nnode), dtype='complex')
    iNR = np.zeros((3,nnode), dtype='complex')
    # Totdal node loads
    sNR = spu*(APQ + AI*np.abs(VNR) + AZ*np.abs(VNR)**2) - 1j*cappu.real + wpu + 1j*vvcpu.real;
    sNR[PH == 0] = 0;
    for ph in range(0,3):
        for k1 in range(0,nnode):
            if np.abs(sNR[ph,k1].real) <= 1e-12:
                sNR[ph,k1] = 0 + sNR[ph,k1].imag
            if np.abs(sNR[ph,k1].imag) <= 1e-12:
                sNR[ph,k1] = sNR[ph,k1].real + 0

    # Total node current
    iNR[PH != 0] = np.conj(sNR[PH != 0]/VNR[PH != 0]); #also needs to be updated...
    iNR[PH == 0] = 0;
    for ph in range(0,3):
        for k1 in range(0,nnode):
            if np.abs(iNR[ph,k1].real) <= 1e-12:
                iNR[ph,k1] = 0 + iNR[ph,k1].imag
            if np.abs(iNR[ph,k1].imag) <= 1e-12:
                iNR[ph,k1] = iNR[ph,k1].real + 0

    print('iNR')
    print(iNR)
    print('sNR')
    print(sNR)

    return VNR, INR, STXNR, SRXNR, iNR, sNR, itercount
