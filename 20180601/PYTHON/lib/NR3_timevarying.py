import numpy as np
from lib.compute_NR3FT_vectorized import compute_NR3FT_vectorized as ft
from lib.compute_NR3JT_vectorized import compute_NR3JT_vectorized as jt
from lib.change_KCL_matrices import change_KCL_matrices
from lib.relevant_openDSS_parameters import relevant_openDSS_parameters
import opendssdirect as dss
import time
import re
def NR3_timevarying(fn, XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, tol, maxiter, der, capacitance, time_delta, H_reg, G_reg):
    dss.run_command('Redirect ' + fn)
    
    vr_no = len(dss.RegControls.AllNames()) #number of voltage regulators
    vr_bus = np.zeros((5, vr_no), dtype = int) 
  
    vr_count = 0 
    vr_lines = 0
    tf_count = 0
    tf_lines = 0

    for vr in range(len(dss.RegControls.AllNames())):
        dss.RegControls.Name(dss.RegControls.AllNames()[vr])  
        for i in range(2): #start / end bus
            dss.Transformers.Name(dss.RegControls.Transformer())     
            bus = dss.CktElement.BusNames()[i].split('.')
            vr_bus[i, vr_count] = int(dss.Circuit.AllBusNames().index(bus[0])) 
        for n in range(len(bus[1:])):                 
            vr_lines += 1                    
            vr_bus[int(bus[1:][n]) + 1, vr_count] = int(bus[1:][n]) # phases that exist                   
        if len(bus) == 1: # unspecified phases, assume all 3 exist
            for n in range(1,4): 
                vr_lines += 1
                vr_bus[n+1, vr_count] = n   
        vr_count += 1 
   
    tf_bus_temp = np.zeros((2, 1))
    for tf in range(len(dss.Transformers.AllNames())):
        dss.Transformers.Name(dss.Transformers.AllNames()[tf]) 
        for i in range(2):     
            bus = dss.CktElement.BusNames()[i].split('.')        
            tf_bus_temp[i] = int(dss.Circuit.AllBusNames().index(bus[0])) 
            # stuff the in and out bus of the tf into an array  
        if not np.size(np.where(vr_bus[0, :] == tf_bus_temp[0])) == 0 and \
        not np.size(np.where(vr_bus[1, :] == tf_bus_temp[1])) == 0:     
            continue #if you have already seen the transformer from regulators, skip
        
        for n in range(len(bus[1:])):                 
            tf_lines += 1                                       
        if len(bus) == 1:
            for _ in range(1,4): #if all phases are assumed to exist
                tf_lines += 1          
        tf_count += 1
   
    nline = len(dss.Lines.AllNames())  
    nnode = len(dss.Circuit.AllBusNames())

    if tol == None:
        tol = 1e-9

    if maxiter == None:
        maxiter = 100

    FT = 1e99
    itercount = 0
    # adjust KCL based on capacitance, DER, and time-varying load
    # if der != 0 or capacitance != 0 or time_delta != -1:  
    #     H, b = change_KCL_matrices(fn, H, g, b, time_delta, der, capacitance)

    # solve power-flow
    t0 = time.time()
    while np.amax(np.abs(FT)) >= 1e-9 and itercount < maxiter:
        print("Iteration number %f" % (itercount))
        FT = ft(XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode, nline, H_reg, G_reg, vr_lines)
        JT = jt(XNR, g_SB, G_KVL, H, g, nnode, nline, H_reg, G_reg, tf_lines, vr_lines)
    
        if JT.shape[0] >= JT.shape[1]:
            XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
        itercount+=1
    t1 = time.time()
    #retrieve relevant parameters for formatting output
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
    print("Time to format VNR")
    t2 = time.time()
    print(t2 - t1)
    print("Time to solve power flow: ")
    print(t1-t0)
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
