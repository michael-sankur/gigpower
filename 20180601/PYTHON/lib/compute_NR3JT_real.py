import numpy as np

def compute_NR3JT_real_function(XNR,network,slackidx,Vslack):

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
    # JT - derivatives of residuals for power flow equations, composed of three parts - see
    # near end of function
    # JTSUBV - derivatives of residuals of slackbus real and imaginary voltage equation
    # components
    # JTKVL - derivatives of residuals of KVL real and imaginary equation components
    # JTKCL - derivatives of residuals of KCL real and imaginary equation components

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

    # Jacobian for slack node voltage
    JSUBV = np.zeros((6,2*3*(nnode + nline)))
    for ph in range(0,3):

        # indexes of the real and imaginary components of the residual
        idxre = 2*ph
        idxim = 2*ph+1

        # indexes of the real and imaginary components voltage
        idxAm = 2*ph*nnode + 2*slackidx
        idxBm = 2*ph*nnode + 2*slackidx+1

        JSUBV[idxre,idxAm] = 1
        JSUBV[idxim,idxBm] = 1

    # Jacobian for KVL across lines (m,n)
    JKVL = np.zeros((2*3*nline,2*3*(nnode + nline)))
    for ph in range(0,3):
        for k1 in range(0,nline):

            # indexes of the real and imaginary components KVL residual
            idxre = 2*ph*nline + 2*k1
            idxim = 2*ph*nline + 2*k1+1

            # if phase does not exist on line
            # I_mn^phi = C_mn^phi + j D_mn^phi = 0
            if LPH[ph,k1] == 0:

                # indexes of C_mn^phi and D_mn^phi
                idxCmn = 2*3*nnode + 2*ph*nline + 2*k1
                idxDmn = 2*3*nnode + 2*ph*nline + 2*k1+1

                # set derivatives of residuals for KVL
                JKVL[idxre,idxCmn] = 1
                JKVL[idxim,idxDmn] = 1

            # if phase does exist on line
            # V_m^phi = V_n^phi + sum_{psi} Z_{mn}^{phi psi} I_{mn}^psi
            # real: A_m^phi = A_n^phi + sum_{psi} r_{mn}^{phi psi} C_{mn}^psi - x_{mn}^{phi psi} D_{mn}^psi
            # imag: B_m^phi = B_n^phi + sum_{psi} r_{mn}^{phi psi} D_{mn}^psi + x_{mn}^{phi psi} C_{mn}^psi
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

                 # set partial derivatives of residuals for KVL

                # real: A_m^phi - A_n^phi - sum_{psi} r_{mn}^{phi,psi} C_{mn}^psi - x_{mn}^{phi,psi} D_{mn}^psi = 0

                # derivatives of real KVl residual with respect to real component of node voltage

                JKVL[idxre,idxAmTx] = 1
                JKVL[idxre,idxAnRx] = -1

                # derivatives of real KVl residual with respect to real component of line current for phases a,b,c
                JKVL[idxre,idxCmna] = -FZpu[ph,3*k1+0].real
                JKVL[idxre,idxCmnb] = -FZpu[ph,3*k1+1].real
                JKVL[idxre,idxCmnc] = -FZpu[ph,3*k1+2].real

                # derivatives of real KVl residual with respect to imag component of line current for phases a,b,c
                JKVL[idxre,idxDmna] = FZpu[ph,3*k1+0].imag
                JKVL[idxre,idxDmnb] = FZpu[ph,3*k1+1].imag
                JKVL[idxre,idxDmnc] = FZpu[ph,3*k1+2].imag

                # imag: B_m^phi - B_n^phi - sum_{psi} r_{mn}^{phi,psi} D_{mn}^psi + x_{mn}^{phi,psi} C_{mn}^psi = 0

                # derivatives of real KVl residual with respect to imag component of node voltage
                JKVL[idxim,idxBmTx] = 1
                JKVL[idxim,idxBnRx] = -1

                # derivatives of imag KVl residual with respect to real component of line current for phases a,b,c
                JKVL[idxim,idxCmna] = -FZpu[ph,3*k1+0].imag
                JKVL[idxim,idxCmnb] = -FZpu[ph,3*k1+1].imag
                JKVL[idxim,idxCmnc] = -FZpu[ph,3*k1+2].imag

                # derivatives of imag KVl residual with respect to imag component of line current for phases a,b,c
                JKVL[idxim,idxDmna] = -FZpu[ph,3*k1+0].real
                JKVL[idxim,idxDmnb] = -FZpu[ph,3*k1+1].real
                JKVL[idxim,idxDmnc] = -FZpu[ph,3*k1+2].real


    # Jacobian for KCL at node m
    JKCL = np.zeros((2*3*(nnode-1),2*3*(nnode + nline)))
    for ph in range(0,3):
        if ph == 0:
            A0 = 1
            B0 = 0
        elif ph == 1:
            A0 = -1/2
            B0 = -np.sqrt(3)/2
        elif ph == 2:
            A0 = -1/2
            B0 = np.sqrt(3)/2
        for k1 in range(1,nnode):
            #if k1 != slackidx:

            # indexes of real and imaginary KCL residual
            idxre = 2*ph*(nnode-1) + 2*(k1-1)
            idxim = 2*ph*(nnode-1) + 2*(k1-1)+1

            # indexes of A_m^phi and B_m^phi
            idxAm = 2*ph*nnode + 2*k1
            idxBm = 2*ph*nnode + 2*k1+1

            # if phase does not exist at node, set V_m^phi = A_m^phi + j B_m^phi = 0
            if NPH[ph,k1] == 0:

                JKCL[idxre,idxAm] = 1
                JKCL[idxim,idxBm] = 1

            # if phase does exist at node
            # sum_{l:(l,m) in Edges} V_m (I_lm^phi)^* = s_m^phi(V_m^phi) + w_m^phi - c_m^phi + sum_{n:(m,n) in Edges} V_m (I_mn^phi)^*
            elif NPH[ph,k1] == 1:

                dA = XNR[idxAm] - A0
                dB = XNR[idxBm] - B0

                gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))])
                gradient_mag_sq = np.array([[2 * A0, 2 * B0]])
                gradient_mag_sq = np.reshape(gradient_mag_sq, (1,2))
                gradient_mag = np.reshape(gradient_mag, (1, 2))

                # # #Using first order Taylor Expansion on magnitude squared (done)
                # #derivates of real KCL residual with respect to real and imag voltage components
                # JKCL[idxre,idxAm] = -spu[ph,k1].real*(AI[ph,k1]* gradient_mag[0] \
                #                     + 2*AZ[ph,k1]*gradient_mag_sq[0][0])
                # JKCL[idxre,idxBm] = -spu[ph,k1].real*(AI[ph,k1] * gradient_mag[1] \
                #                     + 2*AZ[ph,k1]*gradient_mag_sq[0][1])
                #
                # # derivates of imag KCL residual with respect to real and imag voltage components
                # JKCL[idxim,idxAm] = -spu[ph,k1].imag*(AI[ph,k1] * gradient_mag[0] \
                #                                     + 2*AZ[ph,k1]*gradient_mag_sq[0][0])
                # JKCL[idxim,idxBm] = -spu[ph,k1].imag*(AI[ph,k1] * gradient_mag[1] \
                #                                     + 2*AZ[ph,k1]*gradient_mag_sq[0][1])

                #Using second order Taylor Expansion on magnitude squared
                #derivates of real KCL residual with respect to real and imag voltage components
                #hessian_fcn = np.array([[2, 0], [2,0]])
                #
                # dX = np.zeros((2,1))
                # dX[0] = dA
                # dX[1] = dB
                #np.array([dA[0], dB[0]]).T
                # second_order_term = gradient_mag_sq + (np.matmul(hessian_fcn, dX).T)
                #
                # JKCL[idxre,idxAm] = -spu[ph,k1].real*(AI[ph,k1]* gradient_mag[0] \
                #                     + 2*AZ[ph,k1]*second_order_term[0][0])
                # JKCL[idxre,idxBm] = -spu[ph,k1].real*(AI[ph,k1] * gradient_mag[1] \
                #                     + 2*AZ[ph,k1]*second_order_term[0][1])
                #
                # # derivates of imag KCL residual with respect to real and imag voltage components
                # JKCL[idxim,idxAm] = -spu[ph,k1].imag*(AI[ph,k1] * gradient_mag[0] \
                #                                     + 2*AZ[ph,k1]*second_order_term[0][0])
                # JKCL[idxim,idxBm] = -spu[ph,k1].imag*(AI[ph,k1] * gradient_mag[1] \
                #                                     + 2*AZ[ph,k1]*second_order_term[0][1])
                hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
                                    [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])
                #hessian_mag = np.reshape(hessian_mag, (2, 2))
                dX = np.zeros((2,1))
                dX[0] = dA
                dX[1] = dB
                dX = np.reshape(dX, (1, 2))


                second_order_term = (gradient_mag + (dX @ hessian_mag)/2)[0]

                JKCL[idxre,idxAm] = -spu[ph,k1].real*(AI[ph,k1]* second_order_term[0] \
                                    + 2*AZ[ph,k1]*XNR[idxAm])
                JKCL[idxre,idxBm] = -spu[ph,k1].real*(AI[ph,k1] * second_order_term[1] \
                                    + 2*AZ[ph,k1]*XNR[idxBm])

                # derivates of imag KCL residual with respect to real and imag voltage components
                JKCL[idxim,idxAm] = -spu[ph,k1].imag*(AI[ph,k1] * second_order_term[0] \
                                                    + 2*AZ[ph,k1]*XNR[idxAm])
                JKCL[idxim,idxBm] = -spu[ph,k1].imag*(AI[ph,k1] * second_order_term[1] \
                                                    + 2*AZ[ph,k1]*XNR[idxAm])

                # #Using first order Taylor Expansion on magnitude
                # #derivates of real KCL residual with respect to real and imag voltage components
                # JKCL[idxre,idxAm] = -spu[ph,k1].real*(AI[ph,k1]* gradient_mag[0] \
                #                     + 2*AZ[ph,k1]*XNR[idxAm])
                # JKCL[idxre,idxBm] = -spu[ph,k1].real*(AI[ph,k1] * gradient_mag[1] \
                #                     + 2*AZ[ph,k1]*XNR[idxBm])
                #
                # # derivates of imag KCL residual with respect to real and imag voltage components
                # JKCL[idxim,idxAm] = -spu[ph,k1].imag*(AI[ph,k1] * gradient_mag[0] \
                #                                     + 2*AZ[ph,k1]*XNR[idxAm])
                # JKCL[idxim,idxBm] = -spu[ph,k1].imag*(AI[ph,k1] * gradient_mag[1] \
                #                                     + 2*AZ[ph,k1]*XNR[idxBm])


                # #Not Using Taylor Expansion
                # JKCL[idxre,idxAm] = -spu[ph,k1].real*(AI[ph,k1]*XNR[idxAm]*(XNR[idxAm]**2 + XNR[idxBm]**2)**(-1/2) \
                #                                       + 2*AZ[ph,k1]*XNR[idxAm])
                # JKCL[idxre,idxBm] = -spu[ph,k1].real*(AI[ph,k1]*XNR[idxBm]*(XNR[idxAm]**2 + XNR[idxBm]**2)**(-1/2) \
                #                                       + 2*AZ[ph,k1]*XNR[idxBm])
                #
                # # derivates of imag KVL residual with respect to real and imag voltage components
                # JKCL[idxim,idxAm] = -spu[ph,k1].imag*(AI[ph,k1]*XNR[idxAm]*(XNR[idxAm]**2 + XNR[idxBm]**2)**(-1/2) \
                #                                       + 2*AZ[ph,k1]*XNR[idxAm])
                # JKCL[idxim,idxBm] = -spu[ph,k1].imag*(AI[ph,k1]*XNR[idxBm]*(XNR[idxAm]**2 + XNR[idxBm]**2)**(-1/2) \
                #                                       + 2*AZ[ph,k1]*XNR[idxBm])

                # loop through incoming lines to node m - l:(l,m) in Edges
                for k2 in range(0,network.nodes.inlines.shape[0]):

                    # incoming lines connected to node m
                    if inlines[k2,k1] != -1:

                        # indexes of I_lm^phi = C_mn^phi + j D_lm^phi
                        idxClm = 2*3*nnode + 2*ph*nline + 2*inlines[k2,k1]
                        idxDlm = 2*3*nnode + 2*ph*nline + 2*inlines[k2,k1]+1

                        # derivaties of real KCL residual
                        JKCL[idxre,idxAm] = JKCL[idxre,idxAm] + XNR[idxClm]
                        JKCL[idxre,idxBm] = JKCL[idxre,idxBm] + XNR[idxDlm]
                        JKCL[idxre,idxClm] = XNR[idxAm]
                        JKCL[idxre,idxDlm] = XNR[idxBm]

                        # derivaties of imag KCL residual
                        JKCL[idxim,idxAm] = JKCL[idxim,idxAm] - XNR[idxDlm]
                        JKCL[idxim,idxBm] = JKCL[idxim,idxBm] + XNR[idxClm]
                        JKCL[idxim,idxClm] = XNR[idxBm]
                        JKCL[idxim,idxDlm] = -XNR[idxAm]

                # loop through outgoing lines from node m n:(m,n) in Edges
                for k2 in range(0,network.nodes.outlines.shape[0]):

                    # outgoing lines connected to node m
                    if outlines[k2,k1] != -1:

                        # indexes of I_mn^phi = C_mn^phi + j D_mn^phi
                        idxCmn = 2*3*nnode + 2*ph*nline + 2*outlines[k2,k1]
                        idxDmn = 2*3*nnode + 2*ph*nline + 2*outlines[k2,k1]+1

                        # derivaties of real KCL residual
                        JKCL[idxre,idxAm] = JKCL[idxre,idxAm] - XNR[idxCmn]
                        JKCL[idxre,idxBm] = JKCL[idxre,idxBm] - XNR[idxDmn]
                        JKCL[idxre,idxCmn] = -XNR[idxAm]
                        JKCL[idxre,idxDmn] = -XNR[idxBm]

                        # derivaties of imag KCL residual
                        JKCL[idxim,idxAm] = JKCL[idxim,idxAm] + XNR[idxDmn]
                        JKCL[idxim,idxBm] = JKCL[idxim,idxBm] - XNR[idxCmn]
                        JKCL[idxim,idxCmn] = -XNR[idxBm]
                        JKCL[idxim,idxDmn] = XNR[idxAm]

    JT = np.r_[JSUBV, JKVL, JKCL]

    # print('jsubv')
    # print(JSUBV)
    # print('jkvl')
    # print(JKVL)
    # print('jkcl')
    # print(JKCL)
    #
    # a_file = open("not-vectorized.txt", "a+")
    # a_file.write('JSUBV: \n')
    # for row in JSUBV:
    #
    #     np.savetxt(a_file, row)
    # a_file.write('JKVL: \n')
    # for row in JKVL:
    #     np.savetxt(a_file, row)
    # a_file.write('JKCL: \n')
    # for row in JKCL:
    #     np.savetxt(a_file, row)
    #
    # a_file.close()


    return JT
