import numpy as np
import opendssdirect as dss
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

    line_names = dss.Lines.AllNames()

    A_m = np.array([])
    B_m = np.array([])

    C_mn = np.array([])
    D_mn = np.array([])

    R_matrix = np.zeros((nline,9))
    X_matrix = np.zeros((nline,9))

    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])
    kV_base = dss.Bus.kVBase()


    for k1 in range(len(dss.Circuit.AllBusNames())):
        dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k1])
        phases = bus_phase_dict[dss.Circuit.AllBusNames()[k1]]
        volts = dss.Bus.PuVoltage() #get bus1's puVoltage
        a_temp = np.zeros(3)
        b_temp = np.zeros(3)

        count = 0
        for i in range(0, 3):
            if phases[i] == 1:
                a_temp[i] = volts[count]
                b_temp[i] = volts[count+1]
                count = count + 2

        A_m = np.append(A_m, a_temp) #split into re/im parts
        B_m = np.append(B_m, b_temp)


    for k2 in range(len(dss.Lines.AllNames())):
        dss.Lines.Name(dss.Lines.AllNames()[k2]) #set the line

        linecode = dss.Lines.LineCode() #get the linecode
        dss.LineCodes.Name(linecode) #set the linecode

        xmat = dss.LineCodes.Xmatrix() #get the xmat
        rmat = dss.LineCodes.Rmatrix() #get the rmat

        for i in range(len(xmat)):
            X_matrix[k2][i] = xmat[i] #fill x/r where they are shaped like nline x 9 (for 9 components)
        for j in range(len(rmat)):
            R_matrix[k2][i] = rmat[j]

        c_temp = np.zeros(3) #retrieve line current
        d_temp = np.zeros(3)

        for i in range(0, 3): #len(dss.CktElement.Currents()), 2): #get the currents of the line

            c_temp[i] = 0
            d_temp[i] = 0

    #         c_temp[i//2] = np.divide(dss.CktElement.Currents()[i], kV_base) #per unit-ify the currents
    #         d_temp[i//2] = np.divide(dss.CktElement.Currents()[i+1], kV_base)
        C_mn = np.append(C_mn, c_temp)
        D_mn = np.append(D_mn, d_temp)

    X = np.array([]) #make X, should be 2*3*(nline+nnode) long

    for ph in range(0,3):
        for nodes in range(nnode):
            X = np.append(X, A_m[nodes*3 + ph]) #add a, b by node and then phase
            X = np.append(X, B_m[nodes*3 + ph])

    for ph in range(0, 3):
        for lines in range(nline):
            X = np.append(X, C_mn[lines*3 + ph]) #add c, d by line and then phase
            X = np.append(X, D_mn[lines*3 + ph])

    g_SB = np.array([]) #assumes slack bus is at index 0
    sb_idx = [0, 1, 2*nnode, 2*nnode+1, 3*nnode, 3*nnode+1]
    for i in range(len(sb_idx)):
        temp_row = np.zeros(len(X))
        temp_row[sb_idx[i]] = 1
        g_SB = np.append(g_SB, temp_row)
    g_SB = np.reshape(g_SB, (6, len(g_SB) // 6))

    b_SB = np.array([])
    sb_idx = [0, 1, 2*nnode, 2*nnode+1, 3*nnode, 3*nnode+1] #indices of real and im parts of sb voltage
    for i in range(len(sb_idx)):
        b_SB = np.append(b_SB, X[sb_idx[i]])

    FTSUBV = (g_SB @ X) - b_SB

    print('ftsubv \n')
    print(FTSUBV)

    # Residuals for KVL across line (m,n)
    def get_bus_idx(bus):
        k = -1
        for n in range(len(dss.Circuit.AllBusNames())): #iterates over all the buses to see which index corresponds to bus
            if dss.Circuit.AllBusNames()[n] in bus:
                k = n
        return k
    dss.Lines.Name(dss.Lines.AllNames()[0]) #set the line
    bus1 = dss.Lines.Bus1()
    get_bus_idx(bus1)

    def identify_bus_phases(bus): #figures out which phases correspond to the bus
        #returns a list of the r/x matrix places that have those phase/s
        k = np.zeros(3)
        for i in range(1, 4):
            pattern = r"\.%s" % (str(i))
            m = re.findall(pattern, bus)
            if m:
                k[i - 1] = 1
        if np.all(k == np.array([1,0,0])):
            return [0]
        elif np.all(k == np.array([0, 1, 0])):
            return [4]
        elif np.all(k == np.array([0, 0, 1])):
            return [8]
        elif np.all(k == np.array([1, 0, 1])):
            return [0, 6,8]
        elif np.all(k == np.array([1, 1, 0])):
            return [0,3,4]
        elif np.all(k == np.array([0, 1, 1])):
            return [4,7,9]
        else:
            return [0,3,  4, 6, 7, 8]

    G_KVL = np.array([])
    first_template = np.array([1, 0, -1, 0]) #first eqn coeff
    second_template = np.array([0, 1, 0, -1]) #second eqn coeff

    kvl_derivatives = np.zeros((nline , 8))
    kvl_derivatives[:,0] = 1
    kvl_derivatives[:, 1] = -1
    kvl_kvl_derivatives[:, 4] = 1
    derivatives[:, 5] = -1
    for ph in range(0, 3):
        for line in range(len(dss.Lines.AllNames())):
            dss.Lines.Name(dss.Lines.AllNames()[line]) #set the line

            bus1 = dss.Lines.Bus1()
            bus2 = dss.Lines.Bus2()

            b1, b2 = dss.CktElement.BusNames() #the buses on a line should have the same phase
            bus1_phases = identify_bus_phases(b1) #identifies which phase is associated with the bus
            r_count = 0
            x_count = 0

            for e in bus1_phases:
                r_temp = r_count + R_matrix[line, e] #add up all the line resistance comp
                x_temp = x_count + X_matrix[line, e] #add up all the line react comp

            bus1_idx = get_bus_idx(bus1) #get the buses of the line
            bus2_idx = get_bus_idx(bus2)

            temp_row = np.zeros(len(X))

            temp_row[2*nnode*ph + 2*bus1_idx] = 1 #A_m
            temp_row[2*nnode*ph + 2*bus2_idx] = -1 #A_n
            temp_row[2*3*nnode + 2*line*ph] = -r_temp #C_mn
            temp_row[2*3*nnode + 2*line*ph + 1] = x_temp #D_mn
            kvl_derivatives[line, 2] = -r_temp
            kvl_derivatives[line, 3] = x_temp
            G_KVL = np.append(G_KVL, temp_row)

            temp_row = np.zeros(len(X))
            temp_row[2*nnode*ph + 2*bus1_idx + 1] = 1 #B_m
            temp_row[2*nnode*ph + 2*bus2_idx + 1] = -1 #B_n
            temp_row[2*3*nnode + 2*line*ph] = -x_temp #C_mn
            temp_row[2*3*nnode + 2*line*ph + 1] = -r_temp #D_mn
            kvl_derivatives[line, 7] = -r_temp
            kvl_derivatives[line, 6] = -x_temp
            G_KVL = np.append(G_KVL, temp_row)

    G_KVL = np.reshape(G_KVL,(2*3*nline, len(X))) #shape to correct shape

    b_kvl = np.zeros(len(G_KVL))
    FTKVL = np.zeros((2*3*nline,1))
    FTKVL = (G_KVL @ X) - b_kvl
    print('ftkvl \n ')
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
