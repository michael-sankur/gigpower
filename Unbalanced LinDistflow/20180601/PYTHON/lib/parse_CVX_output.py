import numpy as np

def parse_CVX_output_function(X,nvar,network):
    
    # This function parses the output of CVX, which is the variable X

    # INPUT(S)
    # X - optimal solution from CVX
    # nvar - number of variables per phase
    # network - struct containing all pertinent network information,
    # including all other structs
    # base - struct containing base values
    # nodes - struct containing node parameters
    # lines - struct containing line parameters
    # loads - struct containing load parameters
    # caps - struct containing capacitor parameters
    # cons - struct containing controller parameters

    # OUTPUT(S)
    # Eopt - node squared voltage magnitude in 3 x nnode matrix
    # Topt - node voltage angle in 3 x nnode matrix
    # Vopt - node voltage phasor in 3 x nnode matrix
    # Sopt - line complex power in 3 x nline matrix
    # wopt - node DER complex dispatch in 3 x nnode matrix
    # demopt - node complex loads in 3 x nnode matrix
    # sopt - node total complex power in 3 x nnode matrix

    base = network.base
    nodes = network.nodes
    lines = network.lines
    configs = network.configs
    loads = network.loads
    caps = network.caps
    cons = network.cons

    nnode = nodes.nnode
    nline = lines.nline

    # Separate portions of X corresponding to each phase
    Xa = X[0:nvar]
    Xb = X[nvar:2*nvar]
    Xc = X[2*nvar:3*nvar]

    
    # Place in 3 x nvar matrix
    XX = np.r_[Xa.T, Xb.T, Xc.T]
    #print(XX)
    
    XX[np.abs(XX) <= 1e-9] = 0
    XX = np.round(XX, decimals=6)

    # Portion of XX corresponding to the squared voltage magnitude
    Eopt = XX[:,0:nnode]
    #print(Eopt)
    XX = XX[:,nnode:]
    
    # Portion of XX corresponding to the voltage angle
    Topt = XX[:,0:nnode]
    #print(Topt)
    XX = XX[:,nnode:]

    # Compute optimal voltage phasor from Eopt, Topt
    Vopt = np.sqrt(Eopt)*np.exp(1j*np.pi/180*Topt)*nodes.PH
    #print(Vopt)

    Iopt = np.zeros((3,nline), dtype='complex')
    # Compute optimal line current from voltage phasors
    for k1 in range(0,nline):
        #print((Vopt[:,[lines.TXnum[k1]]] - Vopt[:,[lines.RXnum[k1]]]))
        #print(lines.FYpu[:,3*k1:3*k1+3]@(Vopt[:,[lines.TXnum[k1]]] - Vopt[:,[lines.RXnum[k1]]]))
        Iopt[:,[k1]] = lines.FYpu[:,3*k1:3*k1+3]@(Vopt[:,[lines.TXnum[k1]]] - Vopt[:,[lines.RXnum[k1]]])
    Iopt = Iopt*lines.PH

    # Portion of XX corresponding to line real power
    Popt = XX[:,0:nline]
    XX = XX[:,nline:]

    # Portion of XX corresponding to line reactive power
    Qopt = XX[:,0:nline]
    XX = XX[:,nline:]

    # Line complex power
    Sopt = Popt + 1j*Qopt

    # Portion of XX corresponding to node DER real disptach
    uopt = XX[:,0:nnode]
    XX = XX[:,nnode:]

    # Portion of XX corresponding to node DER reactive disptach
    vopt = XX[:,0:nnode]

    # Node complex DER dispatch
    wopt = uopt + 1j*vopt
    #print(wopt)

    # Node loads (demand)
    dopt = loads.spu*(loads.aPQ + loads.aZ*Eopt)*nodes.PH
    #print(demopt)

    # Node total power
    #for k1 = 1:nnode
    #    sopt(:,k1) = sum(Sopt(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
    #        sum(Sopt(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2)
    #end
    sopt = (dopt + wopt - 1j*caps.cappu)*nodes.PH
    
    Eopt = np.round(Eopt, decimals=6)
    Topt = np.round(Topt, decimals=6)
    Vopt = np.round(Vopt, decimals=6)
    Iopt = np.round(Iopt, decimals=6)
    Sopt = np.round(Sopt, decimals=6)
    wopt = np.round(wopt, decimals=6)
    dopt = np.round(dopt, decimals=6)
    sopt = np.round(sopt, decimals=6)
    
    return Eopt, Topt, Vopt, Iopt, Sopt, wopt, dopt, sopt