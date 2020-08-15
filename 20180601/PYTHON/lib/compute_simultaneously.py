import numpy as np
from lib.compute_NR3FT_vectorized import compute_NR3FT_vectorized as ft1
from lib.compute_NR3JT_vectorized import compute_NR3JT_vectorized as ft2
from lib.compute_NR3FT_real import compute_NR3FT_real_function
from lib.compute_NR3JT_real import compute_NR3JT_real_function

from lib.compute_vecmat import compute_vecmat

#from lib.compute_NR3JT_real import compute_NR3JT_real_function
import opendssdirect as dss
import re

def NR3_function(network, fn, slacknode,Vslack,V0,I0,tol=1e-9,maxiter=100):

    dss.run_command('Redirect ' + fn)
    dss.Solution.Solve()
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())
    # This function runs a Newton-Raphson algorithm to solve power flow

    # INPUT(S)
    # network - struct containing all pertinent the network information,
    # including all other structs
    # base - struct containing base values
    # nodes - struct containing node parameters
    # lines - struct containing line parameters
    # loads - struct containing load parameters
    # caps - struct containing capacitor parameters
    # cons - struct containing controller parameters
    # vvc - struct containing vvc paramaters
    # slacknode - index of slack node
    # Vslack - voltage reference for slack node
    # V0 - matrix of node voltages for initializing NR algorthm
    # I0 - matrix of line currents for initializing NR algorthm
    # tol - NR algorithm tolerance
    # maxiter - maximum number of iterations for NR algorthm

    # OUTPUT(S)
    # NRRES - struct containing network states from NR algorthm
    # VNR - Node voltage in 3 x nnode matrix
    # INR - Line current in 3 x nline matrix
    # STXNR - Line power at sending (TX) end in 3 x nline matrix
    # SRXNR - Line power at receiving (RX) end in 3 x nline matrix
    # iNR - Total node current in 3 x nnode matrix - sum of currents from
    # loads, capacitors, and controllers
    # sNR - Total node power in 3 x nnode matrix - sum of the complex power
    # from loads, capacitros and controllers

    # Voltage and current are separated into their real and imaginary parts
    # V_n^phi = A_n^phi + j B_n^phi
    # I_n^phi = C_n^phi + j D_n^phi

    # Voltage and current vectors for a single phase
    # V^phi = [A_1^phi, B_1^phi, A_2^phi, B_2^phi, ... , A_n^phi, B_n^phi]
    # I^phi = [C_1^phi, D_1^phi, C_2^phi, D_2^phi, ... , C_n^phi, D_n^phi]

    # The NR algorithm variable
    # X = [V^a V^b V^c I^a I^b I^c]
    iter = 4 #uncertain how many
    iter_count = 0
    XNR_deep_vec = np.zeros((2*3*(nnode+nline), 1, iter))
    XNR_deep_notvec = np.zeros((2*3*(nnode+nline), 1, iter))
    FT_deep_vec = np.zeros((2*3*(nnode+nline), 1, iter))
    FT_deep_nonvec = np.zeros((2*3*(nnode+nline), 1, iter))
    JT_deep_vec = np.zeros((2*3*(nnode+nline), 2*3*(nnode+nline), iter))
    JT_deep_nonvec = np.zeros((2*3*(nnode+nline), 2*3*(nnode+nline), iter))

    XNR = np.zeros((2*3*(nnode + nline),1))

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
        for k1 in range(0,nnode):
            XNR[(2*3*nnode):] = 0.0*np.ones((6*nline,1))

    elif len(I0) != 0:
        for ph in range(0,3):
            for k1 in range(0,nline):
                XNR[(2*3*nnode) + 2*ph*nline + 2*k1] = I0[ph,k1].real
                XNR[(2*3*nnode) + 2*ph*nline + 2*k1+1] = I0[ph,k1].imag
    if tol == None:
        tol = 1e-9

    if maxiter == None:
        maxiter = 100

    XNR_deep_vec[:, 0, iter_count] = XNR[0]
    XNR_deep_notvec[:, 0, iter_count] = XNR[0]
    iter_count += 1
    FT = 1e99
    itercount = 0
    Vslack = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])


    XNR1, g_SB, b_SB, G_KVL, b_KVL, H, g, b = compute_vecmat(XNR, network, fn, Vslack)
    while np.amax(np.abs(FT)) >= 1e-9 and itercount < maxiter:
        FT1 = ft1(XNR1, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode)
        FT_deep_vec[:, :, iter_count] = FT1
        JT1 = ft2(XNR1, g_SB, G_KVL, H, g, nnode, nline)
        JT_deep_vec[:, :, iter_count] = JT1


        FT = compute_NR3FT_real_function(XNR,network,slacknode,Vslack)
        #print(FT)
        FT_deep_nonvec[:, :, iter_count] = FT
        JT = compute_NR3JT_real_function(XNR,network,slacknode,Vslack)
        JT_deep_nonvec[:, :, iter_count] = JT
        iter_count += 1
        if iter_count == 4:
            return FT_deep_vec, FT_deep_nonvec, JT_deep_vec, JT_deep_nonvec, XNR_deep_vec, XNR_deep_notvec

        if JT.shape[0] >= JT.shape[1]:
            #XNR = XNR - inv(JT.'*JT)*JT.'*FT;
            #print(np.linalg.eigvals(JT))
            #print(np.linalg.eigvals(JT.T@JT))
            XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
        if JT1.shape[0] >= JT1.shape[1]:

            #XNR = XNR - inv(JT.'*JT)*JT.'*FT;
            #print(np.linalg.eigvals(JT))
            #print(np.linalg.eigvals(JT.T@JT))
            XNR1 = XNR1 - np.linalg.inv(JT1.T@JT1)@JT1.T@FT1

        #XNR_deep_notvec[:, 0, iter_count] = np.reshape(XNR1, (len(XNR), 1))
        XNR_deep_vec[:, 0, iter_count] = XNR1[0]
        XNR_deep_notvec[:, 0, iter_count] = XNR[0]

        #print(itercount)
        itercount+=1

        #pass
    #print(XNR)

    #THIS NEEDS TO BE ADJUSTED Because I got rid of the network.

    # remap XNR to VNR, INR, STXNR, SRXNR, iNR, sNR
    # VNR = XNR(1:2:2*3*nnode-1).' + 1j*XNR(2:2:2*3*nnode).';
    # VNR = np.zeros((3,nnode), dtype='complex')
    # for ph in range(0,3):
    #     for k1 in range(0,nnode):
    #         VNR[ph,k1] = XNR[2*ph*nnode + 2*k1] + 1j*XNR[2*ph*nnode + 2*k1+1]
    #         if np.abs(VNR[ph,k1].real) <= 1e-12:
    #             VNR[ph,k1] = 0 + VNR[ph,k1].imag
    #         if np.abs(VNR[ph,k1].imag) <= 1e-12:
    #             VNR[ph,k1] = VNR[ph,k1].real + 0
    # VNR[network.nodes.PH == 0] = 0
    # XNR = XNR[2*3*nnode:]
    #
    #
    # # INR = XNR(2*3*nnode+1:2:2*3*nnode+2*3*nline-1) + 1j*XNR(2*3*nnode+2:2:2*3*nnode+2*3*nline)
    # INR = np.zeros((3,nline), dtype='complex')
    # for ph in range(0,3):
    #     for k1 in range(0,nline):
    #         INR[ph,k1] = XNR[2*ph*nline + 2*k1] + 1j*XNR[2*ph*nline + 2*k1+1]
    #         if np.abs(INR[ph,k1].real) <= 1e-12:
    #             INR[ph,k1] = 0 + INR[ph,k1].imag
    #         if np.abs(INR[ph,k1].imag) <= 1e-12:
    #             INR[ph,k1] = INR[ph,k1].real + 0
    #
    # # STXNR_n^phi = V_m^phi (I_mn^phi)^*
    # # SRXNR_n^phi = V_n^phi (I_mn^phi)^*
    # STXNR = np.zeros((3,nnode), dtype='complex')
    # SRXNR = np.zeros((3,nnode), dtype='complex')
    # for ph in range(0,3):
    #     for k1 in range(0,nline):
    #         STXNR[ph,k1] = VNR[ph,network.lines.TXnum[k1]]*np.conj(INR[ph,k1])
    #         if np.abs(STXNR[ph,k1].real) <= 1e-12:
    #             STXNR[ph,k1] = 0 + STXNR[ph,k1].imag
    #         if np.abs(STXNR[ph,k1].imag) <= 1e-12:
    #             STXNR[ph,k1] = STXNR[ph,k1].real + 0
    #         SRXNR[ph,k1] = VNR[ph,network.lines.RXnum[k1]]*np.conj(INR[ph,k1]) #needs to be updated
    #         if np.abs(SRXNR[ph,k1].real) <= 1e-12:
    #             SRXNR[ph,k1] = 0 + SRXNR[ph,k1].imag
    #         if np.abs(SRXNR[ph,k1].imag) <= 1e-12:
    #             SRXNR[ph,k1] = SRXNR[ph,k1].real + 0
    #
    #
    # # print(STXNR)
    # # print(SRXNR)
    #
    # sNR = np.zeros((3,nnode), dtype='complex')
    # iNR = np.zeros((3,nnode), dtype='complex')
    # # Totdal node loads
    # sNR = spu*(APQ + AI*np.abs(VNR) + AZ*np.abs(VNR)**2) - 1j*cappu.real + wpu + 1j*vvcpu.real;
    # sNR[network.nodes.PH == 0] = 0;
    # for ph in range(0,3):
    #     for k1 in range(0,nnode):
    #         if np.abs(sNR[ph,k1].real) <= 1e-12:
    #             sNR[ph,k1] = 0 + sNR[ph,k1].imag
    #         if np.abs(sNR[ph,k1].imag) <= 1e-12:
    #             sNR[ph,k1] = sNR[ph,k1].real + 0
    # # Total node current
    # iNR[network.nodes.PH != 0] = np.conj(sNR[network.nodes.PH != 0]/VNR[network.nodes.PH != 0]); #also needs to be updated...
    # iNR[network.nodes.PH == 0] = 0;
    # for ph in range(0,3):
    #     for k1 in range(0,nnode):
    #         if np.abs(iNR[ph,k1].real) <= 1e-12:
    #             iNR[ph,k1] = 0 + iNR[ph,k1].imag
    #         if np.abs(iNR[ph,k1].imag) <= 1e-12:
    #             iNR[ph,k1] = iNR[ph,k1].real + 0


    return FT_deep_vec, FT_deep_nonvec, JT_deep_vec, JT_deep_nonvec, XNR_deep_vec, XNR_deep_nonvec
