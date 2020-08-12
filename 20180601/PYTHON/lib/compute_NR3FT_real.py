import numpy as np

def compute_NR3FT_real_function(XNR,network,slackidx,Vslack):

    # Michael Sankur - msankur@lbl.gov
    # 2018.01.01

    # This function computes the residuals of power flow equations

    # INPUT(S)
    # network - struct containing all pertinent the network information,
    # including all other structs
    # base - struct containing base values
    # nodes - struct containing node parameters
    # lines - struct containing line parameters
    # loads - struct containing load parameters
    # caps - struct containing capacitor parameters
    # cons - struct containing controller parameters
    # vvc - struct containing vvc parameters
    # slackidx - index of slack node
    # slackVnom - voltage reference for slack node

    # OUTPUT(S)
    # FT - Residuals for power flow equations, composed of three parts - see
    # near end of function
    # FTSUBV - residuals of slackbus real and imaginary voltage equation
    # components
    # FTKVL - residuals of KVL real and imaginary equation components
    # FTKCL - residuals of KCL real and imaginary equation components

    # slackidx is the node index of the slack bus, which is assigned a fixed
    # voltage reference of slackVnom.

    # Voltage and current are separated into their real and imaginary parts
    # V_n^phi = A_n^phi + j B_n^phi
    # I_n^phi = C_n^phi + j D_n^phi

    # Voltage and current vectors for a single phase
    # V^phi = [A_1^phi, B_1^phi, A_2^phi, B_2^phi, ... , A_n^phi, B_n^phi]
    # I^phi = [C_1^phi, D_1^phi, C_2^phi, D_2^phi, ... , C_n^phi, D_n^phi]

    # The NR algorithm variable
    # X = [V^a V^b V^c I^a I^b I^c]

    '''
    base = network.base
    nodes = network.nodes
    lines = network.lines
    configs = network.configs
    loads = network.loads
    caps = network.caps
    cons = network.cons
    vvc = network.vvc
    '''

    # node parameters
    nnode = network.nodes.nnode
    NPH = network.nodes.PH
    inlines = network.nodes.inlines
    innodes = network.nodes.innodes
    outlines = network.nodes.outlines
    outnodes = network.nodes.outnodes

    # line paramters
    nline = network.lines.nline
    LPH = network.lines.PH
    TXnum = network.lines.TXnum
    RXnum = network.lines.RXnum
    FZpu = network.lines.FZpu
    FRpu = network.lines.FRpu
    FXpu = network.lines.FXpu

    # load parameters
    spu = network.loads.spu
    APQ = network.loads.aPQ
    AI = network.loads.aI
    AZ = network.loads.aZ

    # capacitor paramters
    cappu = network.caps.cappu

    # controller parameters
    wpu = network.cons.wpu

    # vvc parameters
    vvcpu = network.vvc.vvcpu


    # Residuals for slack node voltage
    FTSUBV = np.zeros((6,1))
    FTSUBV[0] = XNR[2*slackidx] - Vslack[0].real
    FTSUBV[1] = XNR[2*slackidx+1] - Vslack[0].imag
    FTSUBV[2] = XNR[2*nnode+2*slackidx] - Vslack[1].real
    FTSUBV[3] = XNR[2*nnode+2*slackidx+1] - Vslack[1].imag
    FTSUBV[4] = XNR[4*nnode+2*slackidx] - Vslack[2].real
    FTSUBV[5] = XNR[4*nnode+2*slackidx+1] - Vslack[2].imag

    # Residuals for KVL across line (m,n)
    FTKVL = np.zeros((2*3*nline,1))
    for ph in range(0,3):
        for k1 in range(0,nline):

            # indexes of real and imag parts of KVL equation for line (m,n)
            idxre = 2*ph*nline + 2*k1
            idxim = 2*ph*nline + 2*k1+1

            # if phase does not exist on line
            # I_mn^phi = C_mn^phi + j D_mn^phi = 0
            if LPH[ph,k1] == 0:

                # indexes of C_mn^phi and D_mn^phi
                idxCmn = 2*3*nnode + 2*ph*nline + 2*k1


                idxDmn = 2*3*nnode + 2*ph*nline + 2*k1+1

                # set residuals for KVL
                FTKVL[idxre] = XNR[idxCmn]
                FTKVL[idxim] = XNR[idxDmn]

            # if phase does exist on line
            # V_m^phi = V_n^phi + sum_{psi} Z_{mn}^{phi psi} I_{mn}^psi
            # Real: A_m^phi = A_n^phi + sum_{psi} r_{mn}^{phi psi} C_{mn}^psi - x_{mn}^{phi psi} D_{mn}^psi
            # Imag: B_m^phi = B_n^phi + sum_{psi} r_{mn}^{phi psi} D_{mn}^psi + x_{mn}^{phi psi} C_{mn}^psi
            elif LPH[ph,k1] == 1:

                # indexes of V_m^phi = A_m^phi +j B_m^phi for TX node
                idxAmTx = 2*ph*nnode + 2*TXnum[k1]
                idxBmTx = 2*ph*nnode + 2*TXnum[k1]+1

                # indexes of V_n^phi = A_n^phi +j B_n^phi for RX node
                idxAnRx = 2*ph*nnode + 2*RXnum[k1]
                idxBnRx = 2*ph*nnode + 2*RXnum[k1]+1

                # indexes of I_mn^a = C_mn^a + j D_mn^a
                idxCmna = 2*3*nnode + 2*k1
                idxDmna = 2*3*nnode + 2*k1+1

                # indexes of I_mn^b = C_mn^b + j D_mn^b
                idxCmnb = 2*3*nnode + 2*nline + 2*k1
                idxDmnb = 2*3*nnode + 2*nline + 2*k1+1

                # indexes of I_mn^c = C_mn^c + j D_mn^c
                idxCmnc = 2*3*nnode + 2*2*nline + 2*k1
                idxDmnc = 2*3*nnode + 2*2*nline + 2*k1+1

                # set residuals for KVL

                # real: A_m^phi - A_n^phi - sum_{psi} [ r_{mn}^{phi,psi} C_{mn}^psi ...\
                # - x_{mn}^{phi,psi} D_{mn}^psi ] = 0
                FTKVL[idxre] = XNR[idxAmTx] - XNR[idxAnRx] \
                    - FZpu[ph,3*k1+0].real*XNR[idxCmna] + FZpu[ph,3*k1+0].imag*XNR[idxDmna] \
                    - FZpu[ph,3*k1+1].real*XNR[idxCmnb] + FZpu[ph,3*k1+1].imag*XNR[idxDmnb] \
                    - FZpu[ph,3*k1+2].real*XNR[idxCmnc] + FZpu[ph,3*k1+2].imag*XNR[idxDmnc]

                # imag: B_m^phi - B_n^phi - sum_{psi} r_{mn}^{phi,psi} D_{mn}^psi ... \
                # + x_{mn}^{phi,psi} C_{mn}^psi = 0
                FTKVL[idxim] = XNR[idxBmTx] - XNR[idxBnRx] \
                    - FZpu[ph,3*k1+0].real*XNR[idxDmna] - FZpu[ph,3*k1+0].imag*XNR[idxCmna] \
                    - FZpu[ph,3*k1+1].real*XNR[idxDmnb] - FZpu[ph,3*k1+1].imag*XNR[idxCmnb] \
                    - FZpu[ph,3*k1+2].real*XNR[idxDmnc] - FZpu[ph,3*k1+2].imag*XNR[idxCmnc]
    print(FTKVL)

    # Residuals for KCL at node m
    # This algorithm assumes that the slack bus has a fixed voltage reference,
    # and its power is "floating" and will be resolved. The slack bus is
    # assumed to be the first node, which respresents the transmission line, or
    # substation if the network configuration is as such - see note below
    FTKCL = np.zeros((2*3*(nnode-1),1))
    for ph in range(0,3):
        if ph == 0:
            A0 = 1
            B0 = 0
        elif ph == 1:
            A0 = -1/2
            B0 = -1 * np.sqrt(3)/2
        elif ph == 2:
            A0 = -1/2
            B0 = np.sqrt(3)/2

        for k1 in range(1,nnode):
            #if k1 != slackidx:

            # indexes of real and imag parts of KCL equation for node m
            idxre = 2*ph*(nnode-1) + 2*(k1-1)
            idxim = 2*ph*(nnode-1) + 2*(k1-1)+1

            # indexes of A_m^phi and B_m^phi
            idxAm = 2*ph*nnode + 2*k1
            idxBm = 2*ph*nnode + 2*k1+1

            #print(slackidx)
            #print(ph, k1, idxre, idxAm)

            # if phase does not exist at node, set V_m^phi = A_m^phi + j B_m^phi = 0
            if NPH[ph,k1] == 0:

                FTKCL[idxre] = XNR[idxAm]
                FTKCL[idxim] = XNR[idxBm]

            # if phase does exist at node
            # sum_{l:(l,m) in Edges} A_m (I_lm^phi)^* = s_m^phi(V_m^phi) + w_m^phi - c_m^phi + sum_{n:(m,n) in Edges} A_m (I_mn^phi)^*
            elif NPH[ph,k1] == 1:

                # initialize residual as zero
                FTKCL[idxre] = 0
                FTKCL[idxim] = 0

                # loop through incoming lines to node m - l:(l,m) in Edges
                for k2 in range(0,network.nodes.inlines.shape[0]):

                    # incoming lines connected to node m
                    if inlines[k2,k1] != -1:

                        # indexes of I_lm^phi = C_lm^phi + j D_lm^phi
                        idxClm = 2*3*nnode + 2*ph*nline + 2*inlines[k2,k1]
                        idxDlm = 2*3*nnode + 2*ph*nline + 2*inlines[k2,k1]+1

                        # sum_{l:(l,m) in Edges} A_m (I_lm^phi)^*
                        # real: A_m^phi C_lm^phi + B_m^phi D_lm^phi
                        # imag: -A_m^phi D_lm^phi + B_m^phi C_lm^phi
                        FTKCL[idxre] = FTKCL[idxre] + XNR[idxAm]*XNR[idxClm] + XNR[idxBm]*XNR[idxDlm]
                        FTKCL[idxim] = FTKCL[idxim] - XNR[idxAm]*XNR[idxDlm] + XNR[idxBm]*XNR[idxClm]

                # s_m^phi(V_m^phi) + w_m^phi - c_m^phi
                # real: p_m^phi (A_{PQ,m}^phi + A_{Z,m}^phi ((A_m^phi)^2 + (B_m^phi)^2))) - u_m^phi
                # imag: q_m^phi (A_{PQ,m}^phi + A_{Z,m}^phi ((A_m^phi)^2 + (B_m^phi)^2))) - v_m^phi + c_m^phi

                dA = XNR[idxAm] - A0
                dB = XNR[idxBm] - B0

                dX = np.array([dA[0], dB[0]])
                dX_t = np.array([dA[0], dB[0]]).T

                gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))])
                gradient_mag_sq = np.array([2 *A0, 2 * B0]) #gradient of magnitude squared

                #Hessian_mag = np.array([[(-1/2) * ],\
                #                        []]) some ratchet chain rule

                # # Applying first order Taylor Expansion to Magnitude Squared (done)
                # FTKCL[idxre] = FTKCL[idxre] \
                #     - spu[ph,k1].real*(APQ[ph,k1] + AI[ph,k1]*
                #     ((A0**2+B0**2)**(1/2) + np.matmul(gradient_mag, dX_t)) \
                #     #(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
                #     + AZ[ph,k1]* \
                #     ((A0**2+B0**2) + np.matmul(gradient_mag_sq, dX_t))) \
                #     - wpu[ph,k1].real
                # FTKCL[idxim] = FTKCL[idxim] \
                #     - spu[ph,k1].imag*(APQ[ph,k1] + AI[ph,k1]* \
                #     ((A0**2+B0**2)**(1/2) + np.matmul(gradient_mag, dX_t)) \
                #     #(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
                #     + AZ[ph,k1]*\
                #     ((A0**2 + B0**2) + np.matmul(gradient_mag_sq, dX_t))) \
                #     + cappu[ph,k1].real - wpu[ph,k1].imag

                # # # Applying second order Taylor Expansion to Magnitude Squared (done)
                # FTKCL[idxre] = FTKCL[idxre] \
                #     - spu[ph,k1].real*(APQ[ph,k1] + AI[ph,k1]*
                #     ((A0**2+B0**2)**(1/2) + np.matmul(gradient_mag, dX_t)) \
                #     #(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
                #     + AZ[ph,k1] * \
                #     ((A0**2+B0**2) + np.matmul(gradient_mag_sq, np.array(dX_t)) + \
                #     (1/2) * np.matmul(dX_t,  2 * dX))) \
                #     - wpu[ph,k1].real
                # FTKCL[idxim] = FTKCL[idxim] \
                #     - spu[ph,k1].imag*(APQ[ph,k1] + AI[ph,k1]* \
                #     ((A0**2+B0**2)**(1/2) + np.matmul(gradient_mag, dX_t)) \
                #     #(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
                #     + AZ[ph,k1]* \
                #     ((A0**2+B0**2) + np.matmul(gradient_mag_sq, dX_t) + \
                #     (1/2) * np.matmul(dX_t, 2 * dX)))  \
                #     + cappu[ph,k1].real - wpu[ph,k1].imag

                # # Applying second order Taylor Expansion to Magnitude Squared and Magnitude
                FTKCL[idxre] = FTKCL[idxre] \
                    - spu[ph,k1].real*(APQ[ph,k1] + AI[ph,k1]*
                    ((A0**2+B0**2)**(1/2) + np.matmul(gradient_mag, dX_t)) \
                    #(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
                    + AZ[ph,k1] * \
                    ((A0**2+B0**2) + np.matmul(gradient_mag_sq, np.array(dX_t)) + \
                    (1/2) * np.matmul(dX_t,  2 * dX))) \
                    - wpu[ph,k1].real
                FTKCL[idxim] = FTKCL[idxim] \
                    - spu[ph,k1].imag*(APQ[ph,k1] + AI[ph,k1]* \
                    ((A0**2+B0**2)**(1/2) + np.matmul(gradient_mag, dX_t)) + \
                    # (1/2) * np.matmul(dX_t, Hessian * dX) #how do you second order expand this
                    #(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
                    + AZ[ph,k1]* \
                    ((A0**2+B0**2) + np.matmul(gradient_mag_sq, dX_t) + \
                    (1/2) * np.matmul(dX_t, 2 * dX)))  \
                    + cappu[ph,k1].real - wpu[ph,k1].imag

                # Applying first order Taylor Expansion to the Magnitude
                # FTKCL[idxre] = FTKCL[idxre] \
                #     - spu[ph,k1].real*(APQ[ph,k1] + AI[ph,k1]*
                #     ((A0**2+B0**2)**(1/2) + np.matmul(gradient_mag, np.array([dA, dB]).T[0])) \
                #     #(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
                #     + AZ[ph,k1]*(XNR[idxAm]**2 + XNR[idxBm]**2)) \
                #     - wpu[ph,k1].real
                # FTKCL[idxim] = FTKCL[idxim] \
                #     - spu[ph,k1].imag*(APQ[ph,k1] + AI[ph,k1]* \
                #     ((A0**2+B0**2)**(1/2) + np.matmul(gradient_mag, np.array([dA, dB]).T[0])) \
                #     #(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
                #     + AZ[ph,k1]*(XNR[idxAm]**2 + XNR[idxBm]**2)) \
                #     + cappu[ph,k1].real - wpu[ph,k1].imag

                # Not using Taylor Expansion
#                 FTKCL[idxre] = FTKCL[idxre] \
#                     - spu[ph,k1].real*(APQ[ph,k1] + AI[ph,k1]*(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
#                     + AZ[ph,k1]*(XNR[idxAm]**2 + XNR[idxBm]**2)) \
#                     - wpu[ph,k1].real
#                 FTKCL[idxim] = FTKCL[idxim] \
#                     - spu[ph,k1].imag*(APQ[ph,k1] + AI[ph,k1]*(XNR[idxAm]**2 + XNR[idxBm]**2)**(1/2) \
#                     + AZ[ph,k1]*(XNR[idxAm]**2 + XNR[idxBm]**2)) \
#                     + cappu[ph,k1].real - wpu[ph,k1].imag

                # loop through outgoing lines from node m - n:(m,n) in Edges
                for k2 in range(0,network.nodes.outlines.shape[0]):

                    # outgoing lines connected to node m
                    if outlines[k2,k1] != -1:

                        # indexes of I_mn^phi = C_mn^phi + j D_mn^phi
                        idxCmn = 2*3*nnode + 2*ph*nline + 2*outlines[k2,k1]
                        idxDmn = 2*3*nnode + 2*ph*nline + 2*outlines[k2,k1]+1

                        # sum_{n:(m,n) in Edges} A_m (I_mn^phi)^*
                        # real: -A_m^phi C_mn^phi - B_m^phi D_mn^phi
                        # imag: A_m^phi D_mn^phi - B_m^phi C_nm^phi
                        FTKCL[idxre] = FTKCL[idxre] - XNR[idxAm]*XNR[idxCmn] - XNR[idxBm]*XNR[idxDmn]
                        FTKCL[idxim] = FTKCL[idxim] + XNR[idxAm]*XNR[idxDmn] - XNR[idxBm]*XNR[idxCmn]

    FT = np.r_[FTSUBV, FTKVL, FTKCL]

    return FT
