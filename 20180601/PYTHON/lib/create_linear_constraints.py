import numpy as np

def create_linear_constraints_function(network,sim,slacknode,Vslack):
    
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
    # slackidx - index of slack node
    # slackVnom - voltage reference for slack node

    # OUTPUT(S)
    # VNR - Node voltage in 3 x nnode matrix
    # INR - Line current in 3 x nline matrix
    # STXNR - Line power at sending (TX) end in 3 x nline matrix
    # SRXNR - Line power at receiving (RX) end in 3 x nline matrix
    # iNR - Total node current in 3 x nnode matrix - sum of currents from
    # loads, capacitors, and controllers
    # sNR - Total node power in 3 x nnode matrix - sum of the complex power
    # from loads, capacitros and controllers

    # Vopt and Iopt are the optimal voltage and current from CVX

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
    inmat = network.nodes.inmat
    outmat = network.nodes.outmat

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
    aPQ = network.loads.aPQ
    aI = network.loads.aI
    aZ = network.loads.aZ

    # capacitor paramters
    cappu = network.caps.cappu

    # controller parameters
    wmaxpu = network.cons.wmaxpu
    
    # simulation parameters
    
    # iteration paramaters
    if hasattr(sim, 'Vkm1'):
        Vkm1 = sim.Vkm1
    else:
        Vkm1 = np.ones((3,nnode), dtype='complex')
        Vkm1[1,:] = np.exp(1j*-120*np.pi/180)
        Vkm1[2,:] = np.exp(1j*120*np.pi/180)

    if hasattr(sim, 'Lkm1'):
        Lkm1 = sim.Lkm1
    else:
        Lkm1 = np.zeros((3,nline))

    if hasattr(sim, 'Hkm1'):
        Hkm1 = sim.Hkm1
    else:
        Hkm1 = np.zeros((3,nline))

    # number of variables per phase
    nvar = nnode + nnode + nline + nline + nnode + nnode
    
    Aeq = np.zeros((1,3*nvar))
    beq = np.zeros((1,1))

    Aineq = np.zeros((1,3*nvar))
    bineq = np.zeros((1,1))
        

    ####################################################################################################
    # Network slack node voltage magnitude constraints
    ####################################################################################################

    # Set network slack node voltage magnitude to specified slack voltage
    
    tempE = np.zeros((1,nnode))
    tempE[0,slacknode] = 1
    tempT = np.zeros((1,nnode))
    tempP = np.zeros((1,nline))
    tempQ = np.zeros((1,nline))
    tempu = np.zeros((1,nnode))
    tempv = np.zeros((1,nnode))

    for ph in range(0,3):
        templine = np.c_[np.zeros((1,ph*nvar)), \
                         tempE, tempT, tempP, tempQ, tempu, tempv, \
                         np.zeros((1,(2-ph)*nvar))]
        Aeq = np.r_[Aeq, templine]

        beq = np.r_[beq, np.abs(Vslack[ph])]
        
    ####################################################################################################
    # Network slack node voltage angle
    ####################################################################################################

    # Set network slack node voltage angles to specified slack voltage

    tempE = np.zeros((1,nnode))
    tempT = np.zeros((1,nnode))
    tempT[0,slacknode] = 1
    tempP = np.zeros((1,nline))
    tempQ = np.zeros((1,nline))
    tempu = np.zeros((1,nnode))
    tempv = np.zeros((1,nnode))

    for ph in range(0,3):
        templine = np.c_[np.zeros((1,ph*nvar)), \
                         tempE, tempT, tempP, tempQ, tempu, tempv, \
                         np.zeros((1,(2-ph)*nvar))]
        Aeq = np.r_[Aeq, templine]

        beq = np.r_[beq, np.angle(Vslack[ph], deg='True')]

    ####################################################################################################
    # Voltage magnitude equality and inequality constraints
    ####################################################################################################
    
    # If no phase present at node - Set voltage magnitude to zero - E_n^phi = 0
    # If phase present at node - Add inequality constraint on squared voltage
    # magnitude - 0.95^2 <= E_n^phi <= 1.05^2

    for ph in range(0,3):
        for k1 in range(0,nnode):
            if NPH[ph,k1] == 0:
                #print(k1, ph)
                tempE = np.zeros((1,nnode))
                tempE[0,k1] = 1
                #print(tempE)
                tempT = np.zeros((1,nnode))
                tempP = np.zeros((1,nline))
                tempQ = np.zeros((1,nline))
                tempu = np.zeros((1,nnode))
                tempv = np.zeros((1,nnode))

                templine = np.c_[np.zeros((1,ph*nvar)), \
                                 tempE, tempT, tempP, tempQ, tempu, tempv, \
                                 np.zeros((1,(2-ph)*nvar))]            
                Aeq = np.r_[Aeq, templine]

                beq = np.r_[beq, [[0]]]

            elif NPH[ph,k1] == 1:
                tempE = np.zeros((2,nnode))
                tempE[:,k1] = [1, -1]
                #print(tempE)
                tempT = np.zeros((2,nnode))
                tempP = np.zeros((2,nline))
                tempQ = np.zeros((2,nline))
                tempu = np.zeros((2,nnode))
                tempv = np.zeros((2,nnode))

                templine = np.c_[np.zeros((2,ph*nvar)), \
                                 tempE, tempT, tempP, tempQ, tempu, tempv, \
                                 np.zeros((2,(2-ph)*nvar))]

                Aineq = np.r_[Aineq, templine]

                bineq = np.r_[bineq, [[1.05**2], [-(0.95**2)]]]

                
    ####################################################################################################
    # Power flow equality constraints
    ####################################################################################################

    # If no phase present at node - Set power entering nodes to zero -
    # P_n^phi = 0, Q_n^phi = 0
    # If phase present at node -
    # P_m^phi - u_m^phi - sum_n ( P_n^phi ) = p_n^phi + sum_n Re{L_mn^phi}
    # Q_m^phi - v_m^phi - sum_n ( Q_n^phi ) = q_n^phi - cap_n^phi + sum_n Im{L_mn^phi}
    # p_n^phi = p_n^phi*(A_PQ,n^phi + A_I,n^phi |V_n^phi| + A_Z,n^phi E_n^phi) + u_n^phi
    # q_n^phi = q_n^phi*(A_PQ,n^phi + A_I,n^phi |V_n^phi| + A_Z,n^phi E_n^phi) - c_n^phi + v_n^phi

    tempE = np.zeros((1,nnode))
    tempT = np.zeros((1,nnode))
    tempu = np.zeros((1,nnode))
    tempv = np.zeros((1,nnode))

    for ph in range(0,3):
        for k1 in range(0,nline):
            if LPH[ph,k1] == 0:
                tempPQ = np.zeros((1,nline))
                tempPQ[0,k1] = 1

                templineP = np.c_[np.zeros((1,ph*nvar)), \
                                 tempE, tempT, tempPQ, np.zeros((1,nline)), tempu, tempv, \
                                 np.zeros((1,(2-ph)*nvar))]
                templineQ = np.c_[np.zeros((1,ph*nvar)), \
                                 tempE, tempT, np.zeros((1,nline)), tempPQ, tempu, tempv, \
                                 np.zeros((1,(2-ph)*nvar))]

                Aeq = np.r_[Aeq, templineP, templineQ]            
                beq = np.r_[beq, np.zeros((2,1))]

    for ph in range(0,3):
        for k1 in range(1,nnode):
            if NPH[ph,k1] == 0:
                for k2 in range(0,inmat.shape[0]):
                    if inmat[k2][k1] != -1:
                        tempE = np.zeros((1,nnode))
                        tempT = np.zeros((1,nnode))
                        tempPQ = np.zeros((1,nline))
                        tempPQ[0,inmat[k2][k1]] = 1
                        tempu = np.zeros((1,nnode))
                        tempv = np.zeros((1,nnode))

                        templineP = np.c_[np.zeros((1,ph*nvar)), \
                                     tempE, tempT, tempPQ, np.zeros((1,nline)), tempu, tempv, \
                                     np.zeros((1,(2-ph)*nvar))]
                        templineQ = np.c_[np.zeros((1,ph*nvar)), \
                                     tempE, tempT, np.zeros((1,nline)), tempPQ, tempu, tempv, \
                                     np.zeros((1,(2-ph)*nvar))]

                        Aeq = np.r_[Aeq, templineP, templineQ]
                        beq = np.r_[beq, np.zeros((2,1))]

                for k2 in range(0,outmat.shape[0]):
                    if outmat[k2][k1] != -1:
                        tempE = np.zeros((1,nnode))
                        tempT = np.zeros((1,nnode))
                        tempPQ = np.zeros((1,nline))
                        tempPQ[0,outmat[k2][k1]] = 1
                        tempu = np.zeros((1,nnode))
                        tempv = np.zeros((1,nnode))

                        templineP = np.c_[np.zeros((1,ph*nvar)), \
                                     tempE, tempT, tempPQ, np.zeros((1,nline)), tempu, tempv, \
                                     np.zeros((1,(2-ph)*nvar))]
                        templineQ = np.c_[np.zeros((1,ph*nvar)), \
                                     tempE, tempT, np.zeros((1,nline)), tempPQ, tempu, tempv, \
                                     np.zeros((1,(2-ph)*nvar))]

                        Aeq = np.r_[Aeq, templineP, templineQ]            
                        beq = np.r_[beq, np.zeros((2,1))]

            elif NPH[ph,k1] == 1:
                tempE = np.zeros((1,nnode), dtype='complex')
                tempE[0,k1] = -spu[ph,k1]*aZ[ph,k1]
                # tempVmag = np.zeros((1,nnode)) tempVmag(1,k1) = -spu(ph,k1)*aI(ph,k1)
                tempT = np.zeros((1,nnode))
                tempPQ = np.zeros((1,nline))
                for k2 in range(0,inmat.shape[0]):
                    if inmat[k2][k1] != -1:
                        tempPQ[0,inmat[k2][k1]] = 1
                for k2 in range(0,inmat.shape[0]):
                    if outmat[k2][k1] != -1:
                        tempPQ[0,outmat[k2][k1]] = -1
                tempuv = np.zeros((1,nnode))
                tempuv[0,k1] = -1

                templineP = np.c_[np.zeros((1,ph*nvar)), \
                             tempE.real, tempT, tempPQ, np.zeros((1,nline)), tempuv, np.zeros((1,nnode)), \
                             np.zeros((1,(2-ph)*nvar))]
                templineQ = np.c_[np.zeros((1,ph*nvar)), \
                             tempE.imag, tempT, np.zeros((1,nline)), tempPQ, np.zeros((1,nnode)), tempuv, \
                             np.zeros((1,(2-ph)*nvar))]

                Aeq = np.r_[Aeq, templineP, templineQ]            
                beq = np.r_[beq, [[aPQ[ph,k1]*spu[ph,k1].real],[aPQ[ph,k1]*spu[ph,k1].imag - cappu[ph,k1]]]]


    ####################################################################################################
    # Calculate gamma terms for each node
    ####################################################################################################
    
    Vrefnode = 'RX'

    Amn = np.zeros((3,3*nnode), dtype='complex')
    Mmn = np.zeros((3,3*nnode))
    Nmn = np.zeros((3,3*nnode))

    for k1 in range(0,nline):

        if Vrefnode == 'TX':
            knode = TXnum[k1]
        if Vrefnode == 'RX':
            knode = RXnum[k1]

        tempAmn = np.zeros((3,3), dtype='complex')
        for ph1 in range(0,3):
            for ph2 in range(0,3):
                if NPH[ph1,knode] == 1 and NPH[ph2,knode] == 1:
                    tempAmn[ph1,ph2] = Vkm1[ph1,knode]/Vkm1[ph2,knode]
        Amn[:,3*k1+0:3*k1+3] = tempAmn
        #print(Amn[:,3*k1+0:3*k1+3])
        Mmn[:,3*k1+0:3*k1+3] = (tempAmn*FZpu[:,3*k1+0:3*k1+3]).real
        #print(Mmn[:,3*k1+0:3*k1+3])
        Nmn[:,3*k1+0:3*k1+3] = (tempAmn*FZpu[:,3*k1+0:3*k1+3]).imag


    ####################################################################################################
    # Voltage magnitude approximation equality constraints
    ####################################################################################################

    # If no phase present at node - Set voltage at phase to zero - E_n^phi = 0
    # If phase present at node - Propagate voltage up feeder (toward head)
    # -E_m + E_n + 2 M_mn P_n + 2 N_mn Q_n = -H_mn

    tempT = np.zeros((1,nnode))
    tempu = np.zeros((1,nnode))
    tempv = np.zeros((1,nnode))

    # Phase A
    for k1 in range(0,nline):

        txnode = TXnum[k1]
        rxnode = RXnum[k1]

        if NPH[0,txnode] == 1 and NPH[0,rxnode] == 1 and LPH[0,k1] == 1:
            tempE = np.zeros((1,nnode))
            tempE[0,[txnode, rxnode]] = [-1, 1]
            tempPa = np.zeros((1,nline))
            tempPa[0,k1] = 2*Mmn[0,3*k1+0]
            tempPb = np.zeros((1,nline))
            tempPb[0,k1] = 2*Mmn[0,3*k1+1]
            tempPc = np.zeros((1,nline))
            tempPc[0,k1] = 2*Mmn[0,3*k1+2]
            tempQa = np.zeros((1,nline))
            tempQa[0,k1] = -2*Nmn[0,3*k1+0]
            tempQb = np.zeros((1,nline))
            tempQb[0,k1] = -2*Nmn[0,3*k1+1]
            tempQc = np.zeros((1,nline))
            tempQc[0,k1] = -2*Nmn[0,3*k1+2]

            templine = np.c_[tempE, tempT, tempPa, tempQa, tempu, tempv, \
                np.zeros((1,nnode)), tempT, tempPb, tempQb, tempu, tempv, \
                np.zeros((1,nnode)), tempT, tempPc, tempQc, tempu, tempv]
            Aeq = np.r_[Aeq, templine]
            beq = np.r_[beq, [[0]]]
            #beq = np.r_[beq, Hkm1[0,k1]]

    # Phase B
    for k1 in range(0,nline):

        txnode = TXnum[k1]
        rxnode = RXnum[k1]

        if NPH[1,txnode] == 1 and NPH[1,rxnode] == 1 and LPH[1,k1] == 1:
            tempE = np.zeros((1,nnode))
            tempE[0,[txnode, rxnode]] = [-1, 1]
            tempPa = np.zeros((1,nline))
            tempPa[0,k1] = 2*Mmn[1,3*k1+0]
            tempPb = np.zeros((1,nline))
            tempPb[0,k1] = 2*Mmn[1,3*k1+1]
            tempPc = np.zeros((1,nline))
            tempPc[0,k1] = 2*Mmn[1,3*k1+2]
            tempQa = np.zeros((1,nline))
            tempQa[0,k1] = -2*Nmn[1,3*k1+0]
            tempQb = np.zeros((1,nline))
            tempQb[0,k1] = -2*Nmn[1,3*k1+1]
            tempQc = np.zeros((1,nline))
            tempQc[0,k1] = -2*Nmn[1,3*k1+2]

            templine = np.c_[np.zeros((1,nnode)), tempT, tempPa, tempQa, tempu, tempv, \
                tempE, tempT, tempPb, tempQb, tempu, tempv, \
                np.zeros((1,nnode)), tempT, tempPc, tempQc, tempu, tempv]
            Aeq = np.r_[Aeq, templine]
            beq = np.r_[beq, [[0]]]
            #beq = np.r_[beq, Hkm1[1,k1]]

    # Phase C
    for k1 in range(0,nline):

        txnode = TXnum[k1]
        rxnode = RXnum[k1]

        if NPH[2,txnode] == 1 and NPH[2,rxnode] == 1 and LPH[2,k1] == 1:
            tempE = np.zeros((1,nnode))
            tempE[0,[txnode, rxnode]] = [-1, 1]
            tempPa = np.zeros((1,nline))
            tempPa[0,k1] = 2*Mmn[2,3*k1+0]
            tempPb = np.zeros((1,nline))
            tempPb[0,k1] = 2*Mmn[2,3*k1+1]
            tempPc = np.zeros((1,nline))
            tempPc[0,k1] = 2*Mmn[2,3*k1+2]
            tempQa = np.zeros((1,nline))
            tempQa[0,k1] = -2*Nmn[2,3*k1+0]
            tempQb = np.zeros((1,nline))
            tempQb[0,k1] = -2*Nmn[2,3*k1+1]
            tempQc = np.zeros((1,nline))
            tempQc[0,k1] = -2*Nmn[2,3*k1+2]

            templine = np.c_[np.zeros((1,nnode)), tempT, tempPa, tempQa, tempu, tempv, \
                np.zeros((1,nnode)), tempT, tempPb, tempQb, tempu, tempv, \
                tempE, tempT, tempPc, tempQc, tempu, tempv]
            Aeq = np.r_[Aeq, templine]
            beq = np.r_[beq, [[0]]]
            #beq = np.r_[beq, Hkm1[2,k1]]

    ####################################################################################################
    # Voltage angle approximation equality constraints
    ####################################################################################################

    # If no phase present at node - No constraint on phase angle
    # If phase present at node - Propagate voltage up feeder (toward head)
    # |V_m||V_n|(-Theta_m + Theta_n) + N_mn P_n - M_mn Q_n = 0
    # Above is not exact equation, read papers to find it

    tempE = np.zeros((1,nnode))
    tempu = np.zeros((1,nnode))
    tempv = np.zeros((1,nnode))

    # Phase A
    for k1 in range(0,nline):

        txnode = TXnum[k1]
        rxnode = RXnum[k1]

        if NPH[0,txnode] == 1 and NPH[0,rxnode] == 1 and LPH[0,k1] == 1:
            tempT = np.zeros((1,nnode))
            tempT[0,[txnode, rxnode]] = [-1, 1]
            tempT = tempT*np.pi/180
            tempPa = np.zeros((1,nline))
            tempPa[0,k1] = -Nmn[0,3*k1+0]
            tempPb = np.zeros((1,nline))
            tempPb[0,k1] = -Nmn[0,3*k1+1]
            tempPc = np.zeros((1,nline))
            tempPc[0,k1] = -Nmn[0,3*k1+2]
            tempQa = np.zeros((1,nline))
            tempQa[0,k1] = -Mmn[0,3*k1+0]
            tempQb = np.zeros((1,nline))
            tempQb[0,k1] = -Mmn[0,3*k1+1]
            tempQc = np.zeros((1,nline))
            tempQc[0,k1] = -Mmn[0,3*k1+2]

            templine = np.c_[tempE, tempT, tempPa, tempQa, tempu, tempv, \
                tempE, np.zeros((1,nnode)), tempPb, tempQb, tempu, tempv, \
                tempE, np.zeros((1,nnode)), tempPc, tempQc, tempu, tempv]
            Aeq = np.r_[Aeq, templine]
            beq = np.r_[beq, [[0]]]


    # Phase B
    for k1 in range(0,nline):

        txnode = TXnum[k1]
        rxnode = RXnum[k1]

        if NPH[1,txnode] == 1 and NPH[1,rxnode] == 1 and LPH[1,k1] == 1:
            tempT = np.zeros((1,nnode))
            tempT[0,[txnode, rxnode]] = [-1, 1]
            tempT = tempT*np.pi/180
            tempPa = np.zeros((1,nline))
            tempPa[0,k1] = -Nmn[1,3*k1+0]
            tempPb = np.zeros((1,nline))
            tempPb[0,k1] = -Nmn[1,3*k1+1]
            tempPc = np.zeros((1,nline))
            tempPc[0,k1] = -Nmn[1,3*k1+2]
            tempQa = np.zeros((1,nline))
            tempQa[0,k1] = -Mmn[1,3*k1+0]
            tempQb = np.zeros((1,nline))
            tempQb[0,k1] = -Mmn[1,3*k1+1]
            tempQc = np.zeros((1,nline))
            tempQc[0,k1] = -Mmn[1,3*k1+2]

            templine = np.c_[tempE, np.zeros((1,nnode)), tempPa, tempQa, tempu, tempv, \
                tempE, tempT, tempPb, tempQb, tempu, tempv, \
                tempE, np.zeros((1,nnode)), tempPc, tempQc, tempu, tempv]
            Aeq = np.r_[Aeq, templine]
            beq = np.r_[beq, [[0]]]


    # Phase C
    for k1 in range(0,nline):

        txnode = TXnum[k1]
        rxnode = RXnum[k1]

        if NPH[2,txnode] == 1 and NPH[2,rxnode] == 1 and LPH[2,k1] == 1:
            tempT = np.zeros((1,nnode))
            tempT[0,[txnode, rxnode]] = [-1, 1]
            tempT = tempT*np.pi/180
            ttempPa = np.zeros((1,nline))
            tempPa[0,k1] = -Nmn[2,3*k1+0]
            tempPb = np.zeros((1,nline))
            tempPb[0,k1] = -Nmn[2,3*k1+1]
            tempPc = np.zeros((1,nline))
            tempPc[0,k1] = -Nmn[2,3*k1+2]
            tempQa = np.zeros((1,nline))
            tempQa[0,k1] = -Mmn[2,3*k1+0]
            tempQb = np.zeros((1,nline))
            tempQb[0,k1] = -Mmn[2,3*k1+1]
            tempQc = np.zeros((1,nline))
            tempQc[0,k1] = -Mmn[2,3*k1+2]

            templine = np.c_[tempE, np.zeros((1,nnode)), tempPa, tempQa, tempu, tempv, \
                tempE, np.zeros((1,nnode)), tempPb, tempQb, tempu, tempv, \
                tempE, tempT, tempPc, tempQc, tempu, tempv]
            Aeq = np.r_[Aeq, templine]
            beq = np.r_[beq, [[0]]]


    ####################################################################################################
    # Zero controller outputs for non control nodes and nonexistent phases
    ####################################################################################################

    # For all phases at all nodes that are non control, set inverter output to
    # zero - u_n^phi = 0, v_n^phi = 0

    tempE = np.zeros((1,nnode))
    tempT = np.zeros((1,nnode))
    tempP = np.zeros((1,nline))
    tempQ = np.zeros((1,nline))

    for ph in range(0,3):
        for k1 in range(0,nnode):
            if NPH[ph,k1] == 0:
                tempuv = np.zeros((1,nnode))
                tempuv[0,k1] = 1

                templineu = np.c_[np.zeros((1,ph*nvar)), \
                                 tempE, tempT, tempP, tempQ, tempuv, np.zeros((1,nnode)), \
                                 np.zeros((1,(2-ph)*nvar))]
                templinev = np.c_[np.zeros((1,ph*nvar)), \
                             tempE, tempT, tempP, tempQ, np.zeros((1,nnode)), tempuv, \
                             np.zeros((1,(2-ph)*nvar))]

                Aeq = np.r_[Aeq, templineu, templinev]            
                beq = np.r_[beq, np.zeros((2,1))]

    ####################################################################################################
    # Inverter control bounds - circle
    ####################################################################################################

    ncirc = 360
    #print(180/np.pi*np.linspace(0,2*np.pi,ncirc+1))

    for ph in range(0,3):
        for k1 in range(0,nnode):
            if NPH[ph,k1] == 0:
                tempE = np.zeros((1,nnode))
                tempT = np.zeros((1,nnode))
                tempP = np.zeros((1,nline))
                tempQ = np.zeros((1,nline))
                tempuv = np.zeros((1,nnode))
                tempuv[0,k1] = 1

                templineu = np.c_[np.zeros((1,ph*nvar)), \
                                 tempE, tempT, tempP, tempQ, tempuv, np.zeros((1,nnode)), \
                                 np.zeros((1,(2-ph)*nvar))]
                templinev = np.c_[np.zeros((1,ph*nvar)), \
                             tempE, tempT, tempP, tempQ, np.zeros((1,nnode)), tempuv, \
                             np.zeros((1,(2-ph)*nvar))]

                Aeq = np.r_[Aeq, templineu, templinev]            
                beq = np.r_[beq, [[0],[0]]]
            elif NPH[ph,k1] == 1 and wmaxpu[ph,k1] == 0:
                tempE = np.zeros((1,nnode))
                tempT = np.zeros((1,nnode))
                tempP = np.zeros((1,nline))
                tempQ = np.zeros((1,nline))
                tempuv = np.zeros((1,nnode))
                tempuv[0,k1] = 1

                templineu = np.c_[np.zeros((1,ph*nvar)), \
                                 tempE, tempT, tempP, tempQ, tempuv, np.zeros((1,nnode)), \
                                 np.zeros((1,(2-ph)*nvar))]
                templinev = np.c_[np.zeros((1,ph*nvar)), \
                             tempE, tempT, tempP, tempQ, np.zeros((1,nnode)), tempuv, \
                             np.zeros((1,(2-ph)*nvar))]

                Aeq = np.r_[Aeq, templineu, templinev]            
                beq = np.r_[beq, [[0],[0]]]
            elif NPH[ph,k1] == 1 and wmaxpu[ph,k1] != 0:
                tempE = np.zeros((ncirc+1,nnode))
                tempT = np.zeros((ncirc+1,nnode))
                tempP = np.zeros((ncirc+1,nline))
                tempQ = np.zeros((ncirc+1,nline))
                tempu = np.zeros((ncirc+1,nnode))
                tempu[:,k1] = np.cos(np.linspace(0,2*np.pi,ncirc+1))
                tempv = np.zeros((ncirc+1,nnode))
                tempv[:,k1] = np.sin(np.linspace(0,2*np.pi,ncirc+1))

                templine = np.c_[np.zeros((ncirc+1,ph*nvar)), \
                                 tempE, tempT, tempP, tempQ, tempu, tempv, \
                                 np.zeros((ncirc+1,(2-ph)*nvar))]

                Aineq = np.r_[Aineq, templine]            
                bineq = np.r_[bineq, wmaxpu[ph,k1]*np.ones((ncirc+1,1))]
                
    Aineq = Aineq[1:,:]
    bineq = bineq[1:]
    Aeq = Aeq[1:,:]
    beq = beq[1:]

    return nvar, Aineq, bineq, Aeq, beq
