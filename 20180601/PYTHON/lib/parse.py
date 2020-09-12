import numpy as np
from lib.compute_NR3FT_vectorized import compute_NR3FT_vectorized as ft1
from lib.compute_NR3JT_vectorized import compute_NR3JT_vectorized as jt1
from lib.change_KCL_matrices import change_KCL_matrices
from lib.relevant_openDSS_parameters import relevant_openDSS_parameters
import opendssdirect as dss
import time
import re
def NR3_timevarying(fn, XNR):
    dss.run_command('Redirect ' + fn)
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())

    t4 = time.time()
    # retrieve relevant parameters for formatting output
    TXnum, RXnum, PH, spu, APQ, AZ, AI, cappu, wpu, vvcpu = \
        relevant_openDSS_parameters(fn, time_delta)
    t5 = time.time()
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


    t6 = time.time()
    print("\n")
    print('change kcl based on time vaying load', t3-t2)
    print('residual calculations/NR3', t4-t3)
    # print('opendss parameters retrival', t5 - t4)
    # print('filling output',t6 - t5)
    # print("\n\n\n")
    return VNR, INR, STXNR, SRXNR, iNR, sNR, itercount
