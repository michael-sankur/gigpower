import numpy as np

# V_n^phi = A_n^phi + j B_n^phi
# I_n^phi = C_n^phi + j D_n^phi

# V^phi = [A_1^phi, B_1^phi, A_2^phi, B_2^phi, ... , A_n^phi, B_n^phi]
# I^phi = [C_1^phi, D_1^phi, C_2^phi, D_2^phi, ... , C_n^phi, D_n^phi]

# X = [V^a V^b V^c I^a I^b I^c]

def computeJTrealfunction(X,feeder,nodes,lines,configs,loads,caps,controllers,slackidx,Vslack):

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

    JSUBV = np.zeros((6,2*3*(nnode + nline)))
    for ph in range(0,3):

        idxre = 2*ph
        idxim = 2*ph+1

        idxAm = 2*ph*nnode + 2*slackidx
        idxBm = 2*ph*nnode + 2*slackidx+1

        JSUBV[idxre,idxAm] = 1
        JSUBV[idxim,idxBm] = 1

    JKVL = np.zeros((2*3*nline,2*3*(nnode + nline)))
    for ph in range(0,3):
        for k1 in range(0,nline):

            idxre = 2*ph*nline + 2*k1
            idxim = 2*ph*nline + 2*k1+1

            if LPH[ph,k1] == 0:

                idxCmn = 2*3*nnode + 2*ph*nline + 2*k1
                idxDmn = 2*3*nnode + 2*ph*nline + 2*k1+1

                JKVL[idxre,idxCmn] = 1
                JKVL[idxim,idxDmn] = 1

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

                JKVL[idxre,idxAmTx] = 1
                JKVL[idxre,idxAmRx] = -1

                JKVL[idxre,idxCmna] = -FZpu[ph,3*k1+0].real
                JKVL[idxre,idxCmnb] = -FZpu[ph,3*k1+1].real
                JKVL[idxre,idxCmnc] = -FZpu[ph,3*k1+2].real

                JKVL[idxre,idxDmna] = FZpu[ph,3*k1+0].imag
                JKVL[idxre,idxDmnb] = FZpu[ph,3*k1+1].imag
                JKVL[idxre,idxDmnc] = FZpu[ph,3*k1+2].imag

                JKVL[idxim,idxBmTx] = 1
                JKVL[idxim,idxBmRx] = -1

                JKVL[idxim,idxCmna] = -FZpu[ph,3*k1+0].imag
                JKVL[idxim,idxCmnb] = -FZpu[ph,3*k1+1].imag
                JKVL[idxim,idxCmnc] = -FZpu[ph,3*k1+2].imag

                JKVL[idxim,idxDmna] = -FZpu[ph,3*k1+0].real
                JKVL[idxim,idxDmnb] = -FZpu[ph,3*k1+1].real
                JKVL[idxim,idxDmnc] = -FZpu[ph,3*k1+2].real


    JKCL = np.zeros((2*3*(nnode-1),2*3*(nnode + nline)))
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

                if NPH[ph,k1] == 0:

                    JKCL[idxre,idxAm] = 1
                    JKCL[idxim,idxBm] = 1

                elif NPH[ph,k1] == 1:

                    JKCL[idxre,idxAm] = -2*spu[ph,k1].real*AZ[ph,k1]*X[idxAm]
                    JKCL[idxre,idxBm] = -2*spu[ph,k1].real*AZ[ph,k1]*X[idxBm]

                    JKCL[idxim,idxAm] = -2*spu[ph,k1].imag*AZ[ph,k1]*X[idxAm]
                    JKCL[idxim,idxBm] = -2*spu[ph,k1].imag*AZ[ph,k1]*X[idxBm]

                    for k2 in range(0,len(inmat[:,k1])):

                        if inmat[k2,k1] != -1:

                            idxClm = 2*3*nnode + 2*ph*nline + 2*inmat[k2,k1]
                            idxDlm = 2*3*nnode + 2*ph*nline + 2*inmat[k2,k1]+1

                            JKCL[idxre,idxAm] = JKCL[idxre,idxAm] + X[idxClm]
                            JKCL[idxre,idxBm] = JKCL[idxre,idxBm] + X[idxDlm]                    
                            JKCL[idxre,idxClm] = X[idxAm]
                            JKCL[idxre,idxDlm] = X[idxBm]


                            JKCL[idxim,idxAm] = JKCL[idxim,idxAm] - X[idxDlm]                    
                            JKCL[idxim,idxBm] = JKCL[idxim,idxBm] + X[idxClm]                 
                            JKCL[idxim,idxClm] = X[idxBm]
                            JKCL[idxim,idxDlm] = -X[idxAm]

                    for k2 in range(0,len(outmat[:,k1])):

                        if outmat[k2,k1] != -1:

                            idxCmn = 2*3*nnode + 2*ph*nline + 2*outmat[k2,k1]
                            idxDmn = 2*3*nnode + 2*ph*nline + 2*outmat[k2,k1]+1

                            JKCL[idxre,idxAm] = JKCL[idxre,idxAm] - X[idxCmn]
                            JKCL[idxre,idxBm] = JKCL[idxre,idxBm] - X[idxDmn]
                            JKCL[idxre,idxCmn] = -X[idxAm]
                            JKCL[idxre,idxDmn] = -X[idxBm]

                            JKCL[idxim,idxAm] = JKCL[idxim,idxAm] + X[idxDmn]
                            JKCL[idxim,idxBm] = JKCL[idxim,idxBm] - X[idxCmn]
                            JKCL[idxim,idxCmn] = -X[idxBm]
                            JKCL[idxim,idxDmn] = X[idxAm]

    JT = np.r_[JSUBV, JKVL, JKCL]

    return JT
                     