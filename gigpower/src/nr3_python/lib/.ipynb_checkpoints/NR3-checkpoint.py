import numpy as np
from lib.computeFTreal import computeFTrealfunction
from lib.computeJTreal import computeJTrealfunction

# V_n^phi = A_n^phi + j B_n^phi
# I_n^phi = C_n^phi + j D_n^phi

# V^phi = [A_1^phi, B_1^phi, A_2^phi, B_2^phi, ... , A_n^phi, B_n^phi]
# I^phi = [C_1^phi, D_1^phi, C_2^phi, D_2^phi, ... , C_n^phi, D_n^phi]

# X = [V^a V^b V^c I^a I^b I^c]

def NR3function(feeder,nodes,lines,configs,loads,caps,controllers,slackidx,Vslack,Vopt,Iopt,maxiter):

    nnode = nodes.nnode
    nline = lines.nline
    
    XNR = np.zeros((2*3*(nnode + nline),1))
    
    # intialize node voltage portion of XNR
    if len(Vopt) == 0:
        for ph in range(0,3):
            for k1 in range(0,nnode):
                XNR[ph*nnode + 2*k1] = Vslack[ph].real
                XNR[ph*nnode + 2*k1+1] = Vslack[ph].imag
    
    elif len(Vopt) != 0:
        for ph in range(0,3):
            for k1 in range(0,nnode):
                XNR[ph*nnode + 2*k1] = Vopt[ph,k1].real
                XNR[ph*nnode + 2*k1+1] = Vopt[ph,k1].imag
        
    
    # intialize line current portion of XNR
    if len(Iopt) == 0:
        for k1 in range(0,nnode):
            XNR[(2*3*nnode):] = 0.1*np.ones((6*nline,1))
            
    elif len(Iopt) != 0:
        for ph in range(0,3):
            for k1 in range(0,nnode):
                XNR[(2*3*nnode) + ph*nnode + 2*k1] = Iopt[ph,k1].real
                XNR[(2*3*nnode) + ph*nnode + 2*k1+1] = Iopt[ph,k1].imag

    
    FT = 1e99
    itercount = 0
    while np.amax(np.abs(FT)) >= 1e-9 and itercount <= maxiter:
        
        FT = computeFTrealfunction(XNR,feeder,nodes,lines,configs,loads,caps,controllers,slackidx,Vslack)
        #print(FT)

        JT = computeJTrealfunction(XNR,feeder,nodes,lines,configs,loads,caps,controllers,slackidx,Vslack)
        #print(JT)
        
        #print(JT.shape)
        
        if JT.shape[0] >= JT.shape[1]:
            #XNR = XNR - inv(JT.'*JT)*JT.'*FT;
            #print(np.linalg.eigvals(JT))
            #print(np.linalg.eigvals(JT.T@JT))
            XNR = XNR - np.linalg.inv(JT.T@JT)@JT.T@FT
        
        #print(itercount)
        itercount+=1
        
        #pass
    #print(XNR)

    VNR = np.zeros((3,nnode), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nnode):
            VNR[ph,k1] = XNR[ph*nnode + 2*k1] + 1j*XNR[ph*nnode + 2*k1+1]
    # print(VNR)
    
    XNR = XNR[2*3*nnode:]
    
    INR = np.zeros((3,nline), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nline):
            INR[ph,k1] = XNR[ph*nnode + 2*k1] + 1j*XNR[ph*nnode + 2*k1+1]
    # print(INR)
    
    STXNR = np.zeros((3,nnode), dtype='complex')
    SRXNR = np.zeros((3,nnode), dtype='complex')
    for k1 in range(0,nline):
        STXNR[:,k1] = VNR[:,lines.TXnum[k1]]*np.conj(INR[:,k1])
        SRXNR[:,k1] = VNR[:,lines.RXnum[k1]]*np.conj(INR[:,k1])
    # print(STXNR)
    # print(SRXNR)
    
    return VNR, INR, STXNR, SRXNR, itercount