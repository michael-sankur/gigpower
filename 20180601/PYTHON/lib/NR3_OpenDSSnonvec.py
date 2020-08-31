import numpy as np
import opendssdirect as dss
from lib.compute_NR3FT_real_OpenDSSnonvec import compute_NR3FT_real_function
from lib.compute_NR3JT_real_OpenDSSnonvec import compute_NR3JT_real_function

def NR3_function(filename,slacknode,Vslack,V0,I0,tol=1e-9,maxiter=100):

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
    dss.run_command('Redirect ' + str(filename))
    dss.Solution.Solve()


    # node parameters
    nnnode = len(dss.Circuit.AllBusNames())
    inlines = network.nodes.inlines
    innodes = network.nodes.innodes
    outlines = network.nodes.outlines
    outnodes = network.nodes.outnodes

    # line paramters
    nline = len(dss.Lines.AllNames())
    TXnum = network.lines.TXnum
    RXnum = network.lines.RXnum
    FZpu = network.lines.FZpu
    FRpu = network.lines.FRpu
    FXpu = network.lines.FXpu

    # load parameters
    #spu = network.loads.spu_nominal
    #spu = network.loads.spu_n
    spu = network.loads.spu
    APQ = 0.85


    AI = 0.05
    AZ = 0.1

    # capacitor paramters
    cappu = dss.Capacitors.kvar()

    # controller parameters
    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k2]) #set bus
    Sk = dss.CktElement.Powers() #retrieve powers
    wpu = Sk[0]

    # vvc parameters
    vvcpu = Sk[1]



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
            if phases[i] == 1: #need to properly assign voltages based on what phases exist
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
            R_matrix[k2][j] = rmat[j]


        c_temp = np.zeros(3) #retrieve line current
        d_temp = np.zeros(3)

        for i in range(0, 3): #len(dss.CktElement.Currents()), 2): #get the currents of the line

            c_temp[i] = 0
            d_temp[i] = 0

    #         c_temp[i//2] = np.divide(dss.CktElement.Currents()[i], kV_base) #per unit-ify the currents
    #         d_temp[i//2] = np.divide(dss.CktElement.Currents()[i+1], kV_base)
        C_mn = np.append(C_mn, c_temp)
        D_mn = np.append(D_mn, d_temp)

    XNR = np.array([]) #make X, X.shape = (2*3*(nline+nnode), 1)

    for ph in range(0,3):
        for nodes in range(nnode):
            XNR = np.append(XNR, A_m[nodes*3 + ph]) #add a, b by node and then phase
            XNR = np.append(XNR, B_m[nodes*3 + ph])

    for ph in range(0, 3):
        for lines in range(nline):
            XNR = np.append(XNR, C_mn[lines*3 + ph]) #add c, d by line and then phase
            XNR = np.append(XNR, D_mn[lines*3 + ph])


    #
    # # intialize node voltage portion of XNR
    # if V0 == None or len(V0) == 0:
    #     for ph in range(0,3):
    #         for k1 in range(0,nnode):
    #             #print(2*ph*nnode + 2*k1)
    #             XNR[2*ph*nnode + 2*k1] = Vslack[ph].real
    #             XNR[2*ph*nnode + 2*k1+1] = Vslack[ph].imag
    #
    # # If initial V is given (usually from CVX)
    # elif len(V0) != 0:
    #     for ph in range(0,3):
    #         for k1 in range(0,nnode):
    #             XNR[2*ph*nnode + 2*k1] = V0[ph,k1].real
    #             XNR[2*ph*nnode + 2*k1+1] = V0[ph,k1].imag
    #
    #
    # # intialize line current portion of XNR
    # if I0 == None or len(I0) == 0:
    #     for k1 in range(0,nnode):
    #         XNR[(2*3*nnode):] = 0.0*np.ones((6*nline,1))
    #
    # elif len(I0) != 0:
    #     for ph in range(0,3):
    #         for k1 in range(0,nline):
    #             XNR[(2*3*nnode) + 2*ph*nline + 2*k1] = I0[ph,k1].real
    #             XNR[(2*3*nnode) + 2*ph*nline + 2*k1+1] = I0[ph,k1].imag

    if tol == None:
        tol = 1e-9

    if maxiter == None:
        maxiter = 100

    FT = 1e99
    itercount = 0
    while np.amax(np.abs(FT)) >= 1e-9 and itercount < maxiter:

        """
        VNR = np.zeros((3,nnode), dtype='complex')
        for ph in range(0,3):
            for k1 in range(0,nnode):
                VNR[ph,k1] = XNR[2*ph*nnode + 2*k1] + 1j*XNR[2*ph*nnode + 2*k1+1]
        # print(VNR)
        """

        FT = compute_NR3FT_real_function(XNR,network,slacknode,Vslack)
        #print(FT)

        JT = compute_NR3JT_real_function(XNR,network,slacknode,Vslack)
        #print(JT)

        if JT.shape[0] >= JT.shape[1]:
            #XNR = XNR - inv(JT.'*JT)*JT.'*FT;
            #print(np.linalg.eigvals(JT))
            #print(np.linalg.eigvals(JT.T@JT))
            XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT

        #print(itercount)
        itercount+=1

        #pass
    #print(XNR)

    # remap XNR to VNR, INR, STXNR, SRXNR, iNR, sNR
    # VNR = XNR(1:2:2*3*nnode-1).' + 1j*XNR(2:2:2*3*nnode).';
    VNR = np.zeros((3,nnode), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nnode):
            VNR[ph,k1] = XNR[2*ph*nnode + 2*k1] + 1j*XNR[2*ph*nnode + 2*k1+1]
            if np.abs(VNR[ph,k1].real) <= 1e-12:
                VNR[ph,k1] = 0 + VNR[ph,k1].imag
            if np.abs(VNR[ph,k1].imag) <= 1e-12:
                VNR[ph,k1] = VNR[ph,k1].real + 0
    VNR[network.nodes.PH == 0] = 0
    XNR = XNR[2*3*nnode:]


    # INR = XNR(2*3*nnode+1:2:2*3*nnode+2*3*nline-1) + 1j*XNR(2*3*nnode+2:2:2*3*nnode+2*3*nline)
    INR = np.zeros((3,nline), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nline):
            INR[ph,k1] = XNR[2*ph*nline + 2*k1] + 1j*XNR[2*ph*nline + 2*k1+1]
            if np.abs(INR[ph,k1].real) <= 1e-12:
                INR[ph,k1] = 0 + INR[ph,k1].imag
            if np.abs(INR[ph,k1].imag) <= 1e-12:
                INR[ph,k1] = INR[ph,k1].real + 0

    # STXNR_n^phi = V_m^phi (I_mn^phi)^*
    # SRXNR_n^phi = V_n^phi (I_mn^phi)^*
    STXNR = np.zeros((3,nnode), dtype='complex')
    SRXNR = np.zeros((3,nnode), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nline):
            STXNR[ph,k1] = VNR[ph,network.lines.TXnum[k1]]*np.conj(INR[ph,k1])
            if np.abs(STXNR[ph,k1].real) <= 1e-12:
                STXNR[ph,k1] = 0 + STXNR[ph,k1].imag
            if np.abs(STXNR[ph,k1].imag) <= 1e-12:
                STXNR[ph,k1] = STXNR[ph,k1].real + 0
            SRXNR[ph,k1] = VNR[ph,network.lines.RXnum[k1]]*np.conj(INR[ph,k1])
            if np.abs(SRXNR[ph,k1].real) <= 1e-12:
                SRXNR[ph,k1] = 0 + SRXNR[ph,k1].imag
            if np.abs(SRXNR[ph,k1].imag) <= 1e-12:
                SRXNR[ph,k1] = SRXNR[ph,k1].real + 0


    # print(STXNR)
    # print(SRXNR)

    sNR = np.zeros((3,nnode), dtype='complex')
    iNR = np.zeros((3,nnode), dtype='complex')
    # Total node loads
    sNR = spu*(APQ + AI*np.abs(VNR) + AZ*np.abs(VNR)**2) - 1j*cappu + wpu + 1j*vvcpu;
    sNR[network.nodes.PH == 0] = 0;
    for ph in range(0,3):
        for k1 in range(0,nnode):
            if np.abs(sNR[ph,k1].real) <= 1e-12:
                sNR[ph,k1] = 0 + sNR[ph,k1].imag
            if np.abs(sNR[ph,k1].imag) <= 1e-12:
                sNR[ph,k1] = sNR[ph,k1].real + 0
    # Total node current
    iNR[network.nodes.PH != 0] = np.conj(sNR[network.nodes.PH != 0]/VNR[network.nodes.PH != 0]);
    iNR[network.nodes.PH == 0] = 0;
    for ph in range(0,3):
        for k1 in range(0,nnode):
            if np.abs(iNR[ph,k1].real) <= 1e-12:
                iNR[ph,k1] = 0 + iNR[ph,k1].imag
            if np.abs(iNR[ph,k1].imag) <= 1e-12:
                iNR[ph,k1] = iNR[ph,k1].real + 0


    return VNR, INR, STXNR, SRXNR, iNR, sNR, itercount
