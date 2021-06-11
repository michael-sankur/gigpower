import numpy as np

def computeFTrealfunction(X,feeder,nodes,lines,configs,loads,caps,controllers,slackidx,Vslack):

    # node parameters
    nnode = nodes.nnode
    NPH = nodes.PH
    inmat = nodes.inmat
    outmat = nodes.outmat

    # line paramters
    nline = lines.nline
    LPH = lines.PH
    TXnum = lines.TXnum
    RXnum = lines.RXnum
    FZpu = lines.FZpu

    # load parameters
    spu = loads.spu
    APQ = loads.aPQ
    AI = loads.aI
    AZ = loads.aZ

    # capacitor paramters
    cappu = caps.cappu

    # controller parameters
    wpu = controllers.wpu

    # Substation fixed voltage
    
    FTSUBV = np.zeros((6,1))
    FTSUBV[0] = X[2*slackidx] - Vslack[0].real
    FTSUBV[1] = X[2*slackidx+1] - Vslack[0].imag
    FTSUBV[2] = X[2*nnode+2*slackidx] - Vslack[1].real
    FTSUBV[3] = X[2*nnode+2*slackidx+1] - Vslack[1].imag
    FTSUBV[4] = X[4*nnode+2*slackidx] - Vslack[2].real
    FTSUBV[5] = X[4*nnode+2*slackidx+1] - Vslack[2].imag
        
    FTKVL = np.zeros((2*3*nline,1))
    for ph in range(0,3):
        for k1 in range(0,nline):
                    
            idxre = 2*ph*nline + 2*k1
            idxim = 2*ph*nline + 2*k1+1
               
            idxCmn = 2*3*nnode + 2*ph*nline + 2*k1
            idxDmn = 2*3*nnode + 2*ph*nline + 2*k1+1

            if LPH[ph,k1] == 0:

                FTKVL[idxre] = X[idxCmn]
                FTKVL[idxim] = X[idxDmn]

            elif LPH[ph,k1] == 1:
            
                idxAmTx = 2*ph*nnode + 2*TXnum[k1]
                idxBmTx = 2*ph*nnode + 2*TXnum[k1]+1

                idxAmRx = 2*ph*nnode + 2*RXnum[k1]            
                idxBmRx = 2*ph*nnode + 2*RXnum[k1]+1

                idxCmna = 2*3*nnode + 2*k1
                idxDmna = 2*3*nnode + 2*k1+1

                idxCmnb = 2*3*nnode + 2*nline + 2*k1
                idxDmnb = 2*3*nnode + 2*nline + 2*k1+1

                idxCmnc = 2*3*nnode + 2*2*nline + 2*k1
                idxDmnc = 2*3*nnode + 2*2*nline + 2*k1+1

                FTKVL[idxre] = X[idxAmTx] - X[idxAmRx] \
                    - FZpu[ph,3*k1+0].real*X[idxCmna] + FZpu[ph,3*k1+0].imag*X[idxDmna] \
                    - FZpu[ph,3*k1+1].real*X[idxCmnb] + FZpu[ph,3*k1+1].imag*X[idxDmnb] \
                    - FZpu[ph,3*k1+2].real*X[idxCmnc] + FZpu[ph,3*k1+2].imag*X[idxDmnc]

                FTKVL[idxim] = X[idxBmTx] - X[idxBmRx] \
                    - FZpu[ph,3*k1+0].real*X[idxDmna] + FZpu[ph,3*k1+0].imag*X[idxCmna] \
                    - FZpu[ph,3*k1+1].real*X[idxDmnb] + FZpu[ph,3*k1+1].imag*X[idxCmnb] \
                    - FZpu[ph,3*k1+2].real*X[idxDmnc] + FZpu[ph,3*k1+2].imag*X[idxCmnc]
    
    
    FTKCL = np.zeros((2*3*(nnode-1),1))
    for ph in range(0,3):
        for k1 in range(0,nnode):
            if k1 != slackidx:
                
                idxre = 2*ph*(nnode-1) + 2*k1
                idxim = 2*ph*(nnode-1) + 2*k1+1
                if k1 > slackidx:
                    idxre -= 2
                    idxim -= 2

                idxAm = 2*ph*nnode + 2*k1
                idxBm = 2*ph*nnode + 2*k1+1
                
                #print(slackidx)
                #print(ph, k1, idxre, idxAm)

                if NPH[ph,k1] == 0:

                    FTKCL[idxre] = X[idxAm]
                    FTKCL[idxim] = X[idxBm]

                elif NPH[ph,k1] == 1:

                    FTKCL[idxre] = 0
                    FTKCL[idxim] = 0

                    for k2 in range(0,len(inmat[:,k1])):

                        if inmat[k2,k1] != -1:

                            idxClm = 2*3*nnode + 2*ph*nline + 2*inmat[k2,k1]
                            idxDlm = 2*3*nnode + 2*ph*nline + 2*inmat[k2,k1]+1

                            FTKCL[idxre] = FTKCL[idxre] + X[idxAm]*X[idxClm] + X[idxBm]*X[idxDlm]                    
                            FTKCL[idxim] = FTKCL[idxim] - X[idxAm]*X[idxDlm] + X[idxBm]*X[idxClm]

                    FTKCL[idxre] = FTKCL[idxre] \
                        - spu[ph,k1].real*(APQ[ph,k1] + AZ[ph,k1]*(X[idxAm]**2 + X[idxBm]**2)) \
                        - wpu[ph,k1].real
                    FTKCL[idxim] = FTKCL[idxim] \
                        - spu[ph,k1].imag*(APQ[ph,k1] + AZ[ph,k1]*(X[idxAm]**2 + X[idxBm]**2)) \
                        - wpu[ph,k1].imag + cappu[ph,k1]

                    for k2 in range(0,len(outmat[:,k1])):

                        if outmat[k2,k1] != -1:

                            idxCmn = 2*3*nnode + 2*ph*nline + 2*outmat[k2,k1]
                            idxDmn = 2*3*nnode + 2*ph*nline + 2*outmat[k2,k1]+1

                            FTKCL[idxre] = FTKCL[idxre] - X[idxAm]*X[idxCmn] - X[idxBm]*X[idxDmn]
                            FTKCL[idxim] = FTKCL[idxim] + X[idxAm]*X[idxDmn] - X[idxBm]*X[idxCmn]
                    

    FT = np.r_[FTSUBV, FTKVL, FTKCL]
    
    return FT
