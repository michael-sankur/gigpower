import numpy as np
import math

printflag = 1

def network_mapper_function(fp, fn):

    fid = open(fp + fn,'r')

    baselinearray = [];
    nodelinearray = [];
    linelinearray = [];
    configlinearray = [];
    loadlinearray = [];
    capacitorlinearray = [];
    controllerlinearray = [];
    vvclinearray = [];

    for templine in fid:
        templine = templine.replace('\n','')
        templine = templine.split(' ')
        if templine[0] == 'base':
            baselinearray.append(templine[1:])
        if templine[0] == 'node':
            nodelinearray.append(templine[1:])
        if templine[0] == 'line':
            linelinearray.append(templine[1:])
        if templine[0] == 'config':
            configlinearray.append(templine[1:])
        if templine[0] == 'load':
            loadlinearray.append(templine[1:])
        if templine[0] == 'capacitor':
            capacitorlinearray.append(templine[1:])
        if templine[0] == 'controller':
            controllerlinearray.append(templine[1:])
        if templine[0] == 'vvc':
            vvclinearray.append(templine[1:])

    '''
    print('BASE\n', baselinearray, '\n')
    print('NODE\n', nodelinearray, '\n')
    print('LINE\n', linelinearray, '\n')
    print('CONFIG\n', configlinearray, '\n')
    print('LOAD\n', loadlinearray, '\n')
    print('CAPACITOR\n', capacitorlinearray, '\n')
    print('CONTROLLER\n', controllerlinearray, '\n')
    '''

    fid.close()


    #########################
    # Classes
    #########################

    class networkclass():
        pass
    network = networkclass()

    class baseclass():
        pass
    base = baseclass()

    class nodesclass():
        pass
    nodes = nodesclass()

    class linesclass():
        pass
    lines = linesclass()

    class configsclass():
        pass
    configs = configsclass()

    class loadsclass():
        pass
    loads = loadsclass()

    class capsclass():
        pass
    caps = capsclass()

    class consclass():
        pass
    cons = consclass()

    class vvcclass():
        pass
    vvc = vvcclass()


    #########################
    # Base Values
    #########################

    # Base values [V], [VAr], [A], [Ohm]

    # Iterate through base value lines
    for k1 in range(0,len(baselinearray)):
        templine = baselinearray[k1]

        for k2 in range(0,len(templine)):

            # Voltage base value [V]
            if templine[k2].split('=')[0] == 'Vbase':
                Vbase = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'unit':
                if templine[k2].split('=')[1] == 'V':
                    pass
                if templine[k2].split('=')[1] == 'kV':
                    Vbase = 1e3*Vbase
                if templine[k2].split('=')[1] == 'MV':
                    Vbase = 1e6*Vbase

            # Voltage base value [VAr]
            if templine[k2].split('=')[0] == 'Sbase':
                Sbase = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'unit':
                if templine[k2].split('=')[1] == 'VAr':
                    pass
                if templine[k2].split('=')[1] == 'kVAr':
                    Sbase = 1e3*Sbase
                if templine[k2].split('=')[1] == 'MVAr':
                    Sbase = 1e6*Sbase

    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase

    base.Vbase = Vbase
    base.Sbase = Sbase
    base.Ibase = Ibase
    base.Zbase = Zbase

    if printflag == 1:
        print('BASE')
        print('Vbase:', base.Vbase)
        print('Sbase:', base.Sbase)
        print('Ibase:', base.Ibase)
        print('Zbase:', base.Zbase)


    #########################
    # Nodes
    #########################

    # List of possible phases for nodes and corresponding matrix
    phlist = ['a','b','c','ab','bc','ac','abc']
    PHmat = np.array([[1, 0, 0], \
                  [0, 1, 0], \
                  [0, 0, 1], \
                  [1, 1, 0], \
                  [0, 1, 1], \
                  [1, 0, 1], \
                  [1, 1, 1]]).T

    # Number of nodes
    nnode = len(nodelinearray)
    nodes.nnode = nnode
    nodes.nodelist = [None]*nnode
    nodes.phases = [None]*nnode
    nodes.PH = np.zeros((3,nnode), dtype='int')

    # Iterate through node in configuration file. Nodes are added in order of their
    # their placement in the configuration file. The slacknode should be the
    # first node in the configuration file (for the time being)
    for k1 in range(0,len(nodelinearray)):
        templine = nodelinearray[k1]

        for k2 in range(0,len(templine)):

            # Node names
            if templine[k2].split('=')[0] == 'nodename':
                nodes.nodelist[k1] = templine[k2].split('=')[1]

            # Node phase list and node phase matrix
            if templine[k2].split('=')[0] == 'phases':
                nodes.phases[k1] = templine[k2].split('=')[1]
                #print(templine[k2].split('=')[1])
                #print(phlist.index(templine[k2].split('=')[1]))
                #print(PHmat[:,phlist.index(templine[k2].split('=')[1])])
                nodes.PH[:,k1] = PHmat[:,phlist.index(templine[k2].split('=')[1])]

    nodes.FM = np.zeros((nnode,nnode))

    if printflag == 1:
        print()
        print('NODES')
        print('nnode:', nodes.nnode)
        print('nodelist:', nodes.nodelist)
        print('phases:', nodes.phases)
        print('PH:\n', nodes.PH)


    #########################
    # Lines
    #########################

    # Number of lines
    nline = len(linelinearray)
    lines.nline = nline
    lines.TXnode = [None]*nline
    lines.RXnode = [None]*nline
    lines.TXnum = np.zeros((len(linelinearray)), dtype='int')
    lines.RXnum = np.zeros((len(linelinearray)), dtype='int')
    lines.phases = [None]*nline
    lines.PH = np.zeros((3,nline), dtype='int')
    lines.config = [None]*nline
    lines.length = np.zeros(len(linelinearray))

    for k1 in range(0,len(linelinearray)):
        templine = linelinearray[k1]

        for k2 in range(0,len(templine)):
            # Line TX (sending) node
            if templine[k2].split('=')[0] == 'TXnode':
                # TX node name placed into string array
                lines.TXnode[k1] = templine[k2].split('=')[1]
                # TX node number (index) placed in array
                lines.TXnum[k1] = int(nodes.nodelist.index(templine[k2].split('=')[1]))

            # Line RX (receiving) node
            if templine[k2].split('=')[0] == 'RXnode':
                # RX node name placed into string array
                lines.RXnode[k1] = templine[k2].split('=')[1]
                # RX node number (index) placed into array
                lines.RXnum[k1] = int(nodes.nodelist.index(templine[k2].split('=')[1]))

            # Line phase list and matrix
            if templine[k2].split('=')[0] == 'phases':
                lines.phases[k1] = templine[k2].split('=')[1]
                lines.PH[:,k1] = PHmat[:,phlist.index(templine[k2].split('=')[1])]

            # Line configuration array
            if templine[k2].split('=')[0] == 'config':
                lines.config[k1] = templine[k2].split('=')[1]

            # Line length [m]
            if templine[k2].split('=')[0] == 'length':
                lines.length[k1] = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'unit':
                if templine[k2].split('=')[1] == 'm':
                    pass
                if templine[k2].split('=')[1] == 'km':
                    lines.length[k1] = 1000*lines.length[k1]
                if templine[k2].split('=')[1] == 'ft' or templine[k2].split('=')[1] == 'feet':
                    lines.length[k1] = 0.3048*lines.length[k1]
                if templine[k2].split('=')[1] == 'mi' or templine[k2].split('=')[1] == 'mile':
                    lines.length[k1] = 1609.34*lines.length[k1]

            nodes.FM[lines.TXnum[k1],lines.RXnum[k1]] = 1
            nodes.FM[lines.RXnum[k1],lines.TXnum[k1]] = -1

    if printflag == 1:
        print()
        print('LINES')
        print('nline:', lines.nline)
        print('TXnode:', lines.TXnode)
        print('RXnode:',lines.RXnode)
        print('TXnum:', lines.TXnum)
        print('RXnum:', lines.RXnum)
        print('phases:', lines.phases)
        print('PH:\n', lines.PH)
        print('config:', lines.config)
        print('length:', lines.length)

    # Find the lines into and out of nodes

    # Compute number of lines entering and leaving each node
    intemp = np.zeros(nnode)
    outtemp = np.zeros(nnode)
    for k1 in range(0,nline):
        intemp[lines.RXnum[k1]]+=1
        outtemp[lines.TXnum[k1]]+=1

    # Matrix of lines coming into nodes, with nodes as columns, and lines as
    # rows. Entries that are -1 are nonexistant (padded)
    nodes.inlines = -np.ones((int(np.amax(intemp)),nnode), dtype='int')
    nodes.innodes = -np.ones((int(np.amax(intemp)),nnode), dtype='int')
    # Matrix of lines leaving nodes, with nodes as columns, and lines as
    # rows. Entries that are -1 are nonexistant (padded)
    nodes.outlines = -np.ones((int(np.amax(outtemp)),nnode), dtype='int')
    nodes.outnodes = -np.ones((int(np.amax(outtemp)),nnode), dtype='int')
    for k1 in range(0,nodes.nnode):
        incount = 0
        outcount = 0
        for k2 in range(0,lines.nline):
            if k1 == lines.RXnum[k2]:
                nodes.inlines[incount,k1] = k2
                nodes.innodes[incount,k1] = lines.TXnum[k2]
                incount+=1
            if k1 == lines.TXnum[k2]:
                nodes.outlines[outcount,k1] = k2
                nodes.outnodes[outcount,k1] = lines.RXnum[k2]
                outcount+=1

    if printflag == 1:
        print()
        print('NODES + LINES')
        print('inlines:\n', nodes.inlines)
        print('innodes:\n', nodes.innodes)
        print('outlines:\n', nodes.outlines)
        print('outnodes:\n', nodes.outnodes)


    #########################
    # Line Configurations
    #########################

    nconf = len(configlinearray)
    configs.nconf = nconf
    configs.conflist = [None]*nconf
    # Impedance matrices for configs [pu/m]
    configs.FZpupl = np.zeros((3,3*nline), dtype='complex')

    # Iterate through line configurations
    for k1 in range(0,len(configlinearray)):
        templine = configlinearray[k1]

        # Impedance (resistance and reactance) values for current line configuration
        tempFZpupl = np.zeros((3,3))
        for k2 in range(0,len(templine)):
            if templine[k2].split('=')[0] == 'config':
                configs.conflist[k1] = (templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'raa':
                raa = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'xaa':
                xaa = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'rab':
                rab = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'xab':
                xab = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'rac':
                rac = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'xac':
                xac = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'rbb':
                rbb = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'xbb':
                xbb = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'rbc':
                rbc = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'xbc':
                xbc = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'rcc':
                rcc = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'xcc':
                xcc = float(templine[k2].split('=')[1])

        tempFZpupl = np.array([[raa + 1j*xaa, rab + 1j*xab, rac + 1j*xac], \
                            [rab + 1j*xab, rbb + 1j*xbb, rbc + 1j*xbc], \
                            [rac + 1j*xac, rbc + 1j*xbc, rcc + 1j*xcc]])

        #print(k1, configs.conflist[k1], tempFZpl)

        # Impedance per unit length multiplying factor for given unit
        for k2 in range(0,len(templine)):
            if templine[k2].split('=')[0] == 'unit':
                if templine[k2].split('=')[1] == 'pu/m':
                    FZmult = 1
                if templine[k2].split('=')[1] == 'pu/km':
                    FZmult = 1e-3
                if templine[k2].split('=')[1] == 'pu/ft' :
                    FZmult = 0.3048
                if templine[k2].split('=')[1] == 'pu/mi':
                    FZmult = 1/1609.34
                if templine[k2].split('=')[1] == 'ohm/m':
                    FZmult = 1/Zbase
                if templine[k2].split('=')[1] == 'ohm/km':
                    FZmult = 1e-3/Zbase
                if templine[k2].split('=')[1] == 'ohm/ft' or templine[k2].split('=')[1] == 'ohm/feet':
                    FZmult = 0.3048/Zbase
                if templine[k2].split('=')[1] == 'ohm/mi' or templine[k2].split('=')[1] == 'ohm/mile':
                    FZmult = 1/1609.34/Zbase

        # 3x3 impedance per unit length matrix for current line config [pu/m]
        configs.FZpupl[:,3*k1:3*(k1+1)] = FZmult*tempFZpupl


    if printflag == 1:
        print()
        print('CONFIGS')
        print('nconf:', configs.nconf)
        print('conflist:', configs.conflist)
        for k1 in range(configs.nconf):
            print(configs.conflist[k1], '- FZpl:\n', configs.FZpupl[:,3*k1:3*(k1+1)])

    # Reconcile Lines and Configurations

    # Impedance matrices for line [pu]
    lines.FZpu = np.zeros((3,3*nline), dtype='complex')
    # Admittance matrices for line [pu]
    lines.FYpu = np.zeros((3,3*nline), dtype='complex')

    # Iterate through lines
    for k1 in range(0,nline):
        # Find configuration for current line
        confnum = configs.conflist.index(lines.config[k1])
        # Multiply 3x3 per unit length impedance matrix [ohm/m] by line length [m]
        # to obtain 3x3 line impdance matrix [pu]
        tempFZpu = lines.length[k1]*configs.FZpupl[:,3*confnum:3*(confnum+1)].reshape((3,3))

        for kph in range(0,3):
            if lines.PH[kph,k1] == 0:
                tempFZpu[:,kph] = 0
                tempFZpu[kph,:] = 0

        lines.FZpu[:,3*k1:3*(k1+1)] = tempFZpu

        tempFYpu = np.zeros((3,3), dtype='complex')
        if lines.phases[k1] == 'a':
            tempFYpu[0,0] = 1/tempFZpu[0,0]
        elif lines.phases[k1] == 'b':
            tempFYpu[1,1] = 1/tempFZpu[1,1]
        elif lines.phases[k1] == 'c':
            tempFYpu[2,2] = 1/tempFZpu[2,2]
        elif lines.phases[k1] == 'ab':
            tempFZpu = tempFZpu[0:2,0:2]
            tempFYpu[0:2,0:2] = np.linalg.inv(tempFZpu)
        elif lines.phases[k1] == 'bc':
            tempFZpu = tempFZpu[1:3,1:3]
            tempFYpu[1:3,1:3] = np.linalg.inv(tempFZpu)
        elif lines.phases[k1] == 'ac':
            tempFZpu = np.array([[tempFZpu[0,0], tempFZpu[0,2]],[tempFZpu[2,0], tempFZpu[2,2]]])
            tempFYpu = np.linalg.inv(tempFZpu)
            tempFYpu = np.insert(tempFYpu,1,[0, 0],0)
            tempFYpu = np.insert(tempFYpu,1,[0, 0, 0],1)
        elif lines.phases[k1] == 'abc':
            tempFYpu = np.linalg.inv(tempFZpu)
        lines.FYpu[:,3*k1:3*(k1+1)] = tempFYpu

    # Resistance, and reactance matrices for all lines [pu]
    lines.FRpu = lines.FZpu.real
    lines.FXpu = lines.FZpu.imag

    # Conductance and susceptance matrices for lines [pu]
    lines.FGpu = lines.FYpu.real
    lines.FBpu = lines.FYpu.imag

    if printflag == 1:
        print()
        print('IMPEDANCE')
        #print('FZpl:', configs.FZpl)
        for k1 in range(lines.nline):
            print(k1, '- FZpu:\n', lines.FZpu[:,3*k1:3*(k1+1)])
        print('ADMITTANCE')
        #print('FZpl:', configs.FZpl)
        for k1 in range(lines.nline):
            print(k1, '- FYpu:\n', lines.FYpu[:,3*k1:3*(k1+1)])


    #########################
    # Loads
    #########################

    loads.conn = []
    # Load parameters
    # Load real component [pu]
    loads.ppu = np.zeros((3,nodes.nnode))
    # Load reactive component [pu]
    loads.qpu = np.zeros((3,nodes.nnode))
    # Load constant power coefficient
    loads.aPQ = np.zeros((3,nodes.nnode))
    # Load constant current coefficient
    loads.aI = np.zeros((3,nodes.nnode))
    # Load constant impedance coefficient
    loads.aZ = np.zeros((3,nodes.nnode))

    # Iterate through loads in configuration
    for k1 in range(0,len(loadlinearray)):
        templine = loadlinearray[k1]

        for k2 in range(0,len(templine)):
            # Load node
            if templine[k2].split('=')[0] == 'nodename':
                knode = nodes.nodelist.index(templine[k2].split('=')[1].upper())

            if templine[k2].split('=')[0] == 'nodenum':
                knode = int(templine[k2].split('=')[1])


            # Load connection - not used
            if templine[k2].split('=')[0] == 'conn':
                pass

            # Load phase(s)
            if templine[k2].split('=')[0] == 'phases':
                if templine[k2].split('=')[1] == 'a':
                    kph = 0
                elif templine[k2].split('=')[1] == 'b':
                    kph = 1
                elif templine[k2].split('=')[1] == 'c':
                    kph = 2

            # Load type - not used
            if templine[k2].split('=')[0] == 'type':
                pass

            # Load constant power coefficient
            if templine[k2].split('=')[0] == 'apq':
                loads.aPQ[kph,knode] = float(templine[k2].split('=')[1])

            # Load constant current coefficient
            if templine[k2].split('=')[0] == 'ai':
                loads.aI[kph,knode] = float(templine[k2].split('=')[1])

            # Load constant impedance coefficient
            if templine[k2].split('=')[0] == 'az':
                loads.aZ[kph,knode] = float(templine[k2].split('=')[1])

            # Load real component [pu]
            if templine[k2].split('=')[0] == 'real':
                loads.ppu[kph,knode] = float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[1] == 'pu':
                    pass
                if templine[k2+1].split('=')[1] == 'W':
                    loads.ppu[kph,knode] = loads.ppu[kph,knode]/Sbase
                elif templine[k2+1].split('=')[1] == 'kW':
                    loads.ppu[kph,knode] = loads.ppu[kph,knode]*1e3/Sbase
                if templine[k2+1].split('=')[1] == 'MW':
                    loads.ppu[kph,knode] = loads.ppu[kph,knode]*1e6/Sbase

            # Load reactive component [pu]
            if templine[k2].split('=')[0] == 'reac':
                loads.qpu[kph,knode] = float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[1] == 'pu':
                    print('apparently is pu')
                    pass
                if templine[k2+1].split('=')[1] == 'VAr':
                    loads.qpu[kph,knode] = loads.qpu[kph,knode]/Sbase
                elif templine[k2+1].split('=')[1] == 'kVAr':
                    print(loads.qpu[kph,knode])
                    loads.qpu[kph,knode] = loads.qpu[kph,knode]*1e3/Sbase
                if templine[k2+1].split('=')[1] == 'MVAr':
                    loads.qpu[kph,knode] = loads.qpu[kph,knode]*1e6/Sbase

    # Complex loads [pu]
    loads.spu = loads.ppu + 1j*loads.qpu
    loads.spu_nominal = loads.ppu +1j*loads.qpu



#     loads.ppu_n = np.abs(loads.ppu * np.random.normal(0, 1, loads.ppu.shape))
#     loads.qpu_n = np.abs(loads.qpu * np.random.normal(0, 1, loads.ppu.shape))
#     loads.spu_n = loads.ppu_n + loads.qpu_n


    if printflag == 1:
        print()
        print('LOADS')
        print('aPQ:\n', loads.aPQ)
        print('aI:\n', loads.aI)
        print('aZ:\n', loads.aZ)
        print('ppu:\n', loads.ppu)
        print('qpu:\n', loads.qpu)
        print('spu:\n', loads.spu)


    #########################
    # Capacitors
    #########################

    # Capacitor parameters
    # Capacitor matrix [pu]
    caps.cappu = np.zeros((3,nodes.nnode))

    # Iterate through capacitors in configuration file
    for k1 in range(0,len(capacitorlinearray)):
        templine = capacitorlinearray[k1]

        for k2 in range(0,len(templine)):

            # Capacitor node
            if templine[k2].split('=')[0] == 'nodename':
                knode = nodes.nodelist.index(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'nodenum':
                knode = int(templine[k2].split('=')[1])

            # Capacitor connection - not used
            if templine[k2].split('=')[0] == 'conn':
                if templine[k2].split('=')[1] == 'wye':
                    pass
                if templine[k2].split('=')[1] == 'delta':
                    pass

            # Capacitor phase(s)
            if templine[k2].split('=')[0] == 'phase':
                if templine[k2].split('=')[1] == 'a':
                    kph = 0
                elif templine[k2].split('=')[1] == 'b':
                    kph = 1
                elif templine[k2].split('=')[1] == 'c':
                    kph = 2

            # Capacitance [pu]
            if templine[k2].split('=')[0] == 'reac':
                caps.cappu[kph,knode] = float(templine[k2].split('=')[1])
                #caps.cap[kph,knode] = 1j*float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[0] == 'unit':
                    if templine[k2+1].split('=')[1] == 'pu':
                        pass
                    if templine[k2+1].split('=')[1] == 'VAr':
                        caps.cappu[kph,knode] = caps.cappu[kph,knode]/Sbase
                    if templine[k2+1].split('=')[1] == 'kVAr':
                        caps.cappu[kph,knode] = caps.cappu[kph,knode]*1e3/Sbase
                    if templine[k2+1].split('=')[1] == 'MVAr':
                        caps.cappu[kph,knode] = caps.cappu[kph,knode]*1e6/Sbase

    if printflag == 1:
        print()
        print('CAPS')
        print('cappu:\n', caps.cappu)


    #########################
    # Controllers
    #########################

    # Controller parameters
    # Controller apparent power capacity [pu]
    cons.wmaxpu = np.zeros((3,nodes.nnode))

    # ES controller frequency
    cons.fes = np.zeros((3,nodes.nnode))
    # ES controller high pass filter (hpf) frequency
    cons.hpfes = np.zeros((3,nodes.nnode))
    # ES controller low pass filter (lpf) frequency
    cons.lpfes = np.zeros((3,nodes.nnode))
    # ES controller integrator gain
    cons.kintes = np.zeros((3,nodes.nnode))

    for k1 in range(0,len(controllerlinearray)):
        templine = controllerlinearray[k1]

        for k2 in range(0,len(templine)):

            # Controller node
            if templine[k2].split('=')[0] == 'nodename':
                knode = nodes.nodelist.index(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'nodenum':
                knode = int(templine[k2].split('=')[1])

            # Controller connection - not used
            if templine[k2].split('=')[0] == 'conn':
                if templine[k2].split('=')[1] == 'wye':
                    pass
                if templine[k2].split('=')[1] == 'delta':
                    pass

            # Controller phases(s)
            if templine[k2].split('=')[0] == 'phases':
                if templine[k2].split('=')[1] == 'a':
                    kph = 0
                elif templine[k2].split('=')[1] == 'b':
                    kph = 1
                elif templine[k2].split('=')[1] == 'c':
                    kph = 2

            # Controller apparent power capacity [pu]
            if templine[k2].split('=')[0] == 'mag':
                cons.wmaxpu[kph,knode] = float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[0] == 'unit':
                    if templine[k2+1].split('=')[1] == 'pu':
                        pass
                    if templine[k2+1].split('=')[1] == 'VAr':
                        cons.wmaxpu[kph,knode] = cons.wmaxpu[kph,knode]/Sbase
                    if templine[k2+1].split('=')[1] == 'kVAr':
                        cons.wmaxpu[kph,knode] = cons.wmaxpu[kph,knode]*1e3/Sbase
                    if templine[k2+1].split('=')[1] == 'MVAr':
                        cons.wmaxpu[kph,knode] = cons.wmaxpu[kph,knode]*1e6/Sbase

            # Controller ES frequency
            if templine[k2].split('=')[0] == 'fes':
                cons.fes[kph,knode] = float(templine[k2].split('=')[1])

            # Controller ES high pass filter (hpf) frequency
            if templine[k2].split('=')[0] == 'hpfes':
                cons.hpfes[kph,knode] = float(templine[k2].split('=')[1])

            # Controller ES low pass filter (lpf) frequency
            if templine[k2].split('=')[0] == 'lpfes':
                cons.lpfes[kph,knode] = float(templine[k2].split('=')[1])

            # Controller ES integrator gain
            if templine[k2].split('=')[0] == 'kintes':
                cons.kintes[kph,knode] = float(templine[k2].split('=')[1])

    # Controller dispatch [pu]
    cons.wpu = np.zeros((3,nodes.nnode))

    if printflag == 1:
        print()
        print('CONS')
        print('wmaxpu:\n', cons.wmaxpu)
        print('wpu:\n', cons.wpu)
        print('fes:\n', cons.fes)
        print('hpfes:\n', cons.hpfes)
        print('lpfes:\n', cons.lpfes)
        print('kintes:\n', cons.kintes)


    #########################
    # VVC
    #########################


    vvc.state =np.zeros((3,nnode))
    vvc.type =np.zeros((3,nnode))
    # VVC minimum voltage [pu]
    vvc.Vminpu =np.zeros((3,nnode))
    # VVC maximum voltage [pu]
    vvc.Vmaxpu =np.zeros((3,nnode))
    # VVC maximum reactive power [pu]
    vvc.qminpu =np.zeros((3,nnode))
    # VVC maximum reactive power [pu]
    vvc.qmaxpu =np.zeros((3,nnode))

    for k1 in range(0,len(controllerlinearray)):
        templine = controllerlinearray[k1]

        for k2 in range(0,len(templine)):

            # VVC node
            if templine[k2].split('=')[0] == 'nodename':
                knode = nodes.nodelist.index(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'nodenum':
                knode = int(templine[k2].split('=')[1])

            # VVC connection - not used
            if templine[k2].split('=')[0] == 'conn':
                if templine[k2].split('=')[1] == 'wye':
                    pass
                if templine[k2].split('=')[1] == 'delta':
                    pass

            # VVC phases(s)
            if templine[k2].split('=')[0] == 'phases':
                if templine[k2].split('=')[1] == 'a':
                    kph = 0
                elif templine[k2].split('=')[1] == 'b':
                    kph = 1
                elif templine[k2].split('=')[1] == 'c':
                    kph = 2

            # VVC state
            if templine[k2].split('=')[0] == 'state':
                vvc.state[kph,knode] = int(templine[k2].split('=')[1])

            # VVC type
            if templine[k2].split('=')[0] == 'type':
                vvc.type[kph,knode] = int(templine[k2].split('=')[1])

            # VVC minimum voltage [pu]
            if templine[k2].split('=')[0] == 'Vmin':
                vvc.Vminpu[kph,knode] = float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[0] == 'unit':
                    if templine[k2+1].split('=')[1] == 'pu':
                        pass
                    if templine[k2+1].split('=')[1] == 'V':
                        vvc.Vminpu[kph,knode] = vvc.Vminpu[kph,knode]/Vbase
                    if templine[k2+1].split('=')[1] == 'kV':
                        vvc.Vminpu[kph,knode] = vvc.Vminpu[kph,knode]*1e3/Vbase
                    if templine[k2+1].split('=')[1] == 'MV':
                        vvc.Vminpu[kph,knode] = vvc.Vminpu[kph,knode]*1e6/Vbase

            # VVC maximum voltage [pu]
            if templine[k2].split('=')[0] == 'Vmax':
                vvc.Vmaxpu[kph,knode] = float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[0] == 'unit':
                    if templine[k2+1].split('=')[1] == 'pu':
                        pass
                    if templine[k2+1].split('=')[1] == 'V':
                        vvc.Vmaxpu[kph,knode] = vvc.Vmaxpu[kph,knode]/Vbase
                    if templine[k2+1].split('=')[1] == 'kV':
                        vvc.Vmaxpu[kph,knode] = vvc.Vmaxpu[kph,knode]*1e3/Vbase
                    if templine[k2+1].split('=')[1] == 'MV':
                        vvc.Vmaxpu[kph,knode] = vvc.Vmaxpu[kph,knode]*1e6/Vbase

            # VVC minimum reactive power
            if templine[k2].split('=')[0] == 'qmin':
                vvc.qminpu[kph,knode] = float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[0] == 'unit':
                    if templine[k2+1].split('=')[1] == 'pu':
                        pass
                    if templine[k2+1].split('=')[1] == 'VAr':
                        vvc.qminpu[kph,knode] = vvc.qminpu[kph,knode]/Sbase
                    if templine[k2+1].split('=')[1] == 'kVAr':
                        vvc.qminpu[kph,knode] = vvc.qminpu[kph,knode]*1e3/Sbase
                    if templine[k2+1].split('=')[1] == 'MVAr':
                        vvc.qminpu[kph,knode] = vvc.qminpu[kph,knode]*1e6/Sbase


            # VVC maximum reactive power
            if templine[k2].split('=')[0] == 'qmin':
                vvc.qmaxpu[kph,knode] = float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[0] == 'unit':
                    if templine[k2+1].split('=')[1] == 'pu':
                        pass
                    if templine[k2+1].split('=')[1] == 'VAr':
                        vvc.qmaxpu[kph,knode] = vvc.qmaxpu[kph,knode]/Sbase
                    if templine[k2+1].split('=')[1] == 'kVAr':
                        vvc.qmaxpu[kph,knode] = vvc.qmaxpu[kph,knode]*1e3/Sbase
                    if templine[k2+1].split('=')[1] == 'MVAr':
                        vvc.qmaxpu[kph,knode] = vvc.qmaxpu[kph,knode]*1e6/Sbase


    # VVC dispatch [pu]
    vvc.vvcpu = np.zeros((3,nnode))

    if printflag == 1:
        print()
        print('VVC')
        print('state:\n', vvc.state)
        print('type:\n', vvc.type)
        print('Vminpu:\n', vvc.Vminpu)
        print('Vmaxpu:\n', vvc.Vmaxpu)
        print('qminpu:\n', vvc.qminpu)
        print('qmaxpu:\n', vvc.qmaxpu)


    #########################
    # NETWORK
    #########################

    network.base = base
    network.nodes = nodes
    network.lines = lines
    network.configs = configs
    network.loads = loads
    network.caps = caps
    network.cons = cons
    network.vvc = vvc

    return network
