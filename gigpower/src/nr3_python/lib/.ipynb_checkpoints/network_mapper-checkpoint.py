import numpy as np
import math

phlist = ['a','b','c','ab','bc','ac','abc']
PHmat = np.array([[1, 0, 0], \
                  [0, 1, 0], \
                  [0, 0, 1], \
                  [1, 1, 0], \
                  [0, 1, 1], \
                  [1, 0, 1], \
                  [1, 1, 1]]).T

def network_mapper_function(fp, fn):

    fid = open(fp + fn,'r')

    baselinearray = [];
    nodelinearray = [];
    linelinearray = [];
    configlinearray = [];
    loadlinearray = [];
    capacitorlinearray = [];
    controllerlinearray = [];

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

    # Base Values

    class networkclass():
        pass
    network = networkclass()

    for k1 in range(0,len(baselinearray)):    
        templine = baselinearray[k1]

        for k2 in range(0,len(templine)):
            #print(templine[k2])
            if templine[k2].split('=')[0] == 'Vbase':
                network.Vbase = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'unit':
                if templine[k2].split('=')[1] == 'V':
                    pass
                if templine[k2].split('=')[1] == 'kV':
                    network.Vbase = 1e3*network.Vbase
                if templine[k2].split('=')[1] == 'MV':
                    network.Vbase = 1e6*network.Vbase

            if templine[k2].split('=')[0] == 'Sbase':
                network.Sbase = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'unit':
                if templine[k2].split('=')[1] == 'VAr':
                    pass
                if templine[k2].split('=')[1] == 'kVAr':
                    network.Sbase = 1e3*network.Sbase
                if templine[k2].split('=')[1] == 'MVAr':
                    network.Sbase = 1e6*network.Sbase

    network.Ibase = network.Sbase/network.Vbase
    network.Zbase = network.Vbase/network.Ibase

    print(network.Vbase, network.Sbase, network.Ibase, network.Zbase)


    # Nodes

    class nodesclass():
        pass
    nodes = nodesclass()

    nnode = len(nodelinearray)
    nodes.nnode = nnode
    nodes.nodelist = [None]*nnode
    nodes.phases = [None]*nnode
    nodes.PH = np.zeros((3,nnode), dtype='int')

    for k1 in range(0,len(nodelinearray)):
        templine = nodelinearray[k1]

        for k2 in range(0,len(templine)):
            #print(templine[k2])
            if templine[k2].split('=')[0] == 'nodename':
                nodes.nodelist[k1] = templine[k2].split('=')[1]
            if templine[k2].split('=')[0] == 'phases':
                nodes.phases[k1] = templine[k2].split('=')[1]
                #print(templine[k2].split('=')[1])
                #print(phlist.index(templine[k2].split('=')[1]))
                #print(PHmat[:,phlist.index(templine[k2].split('=')[1])])
                nodes.PH[:,k1] = PHmat[:,phlist.index(templine[k2].split('=')[1])]

    print(nodes.nodelist)
    print(nodes.phases)
    print(nodes.PH)

    # Lines

    class linesclass():
        pass
    lines = linesclass()

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
            if templine[k2].split('=')[0] == 'TXnode':
                lines.TXnode[k1] = templine[k2].split('=')[1]
                lines.TXnum[k1] = int(nodes.nodelist.index(templine[k2].split('=')[1]))
            if templine[k2].split('=')[0] == 'RXnode':
                lines.RXnode[k1] = templine[k2].split('=')[1]
                lines.RXnum[k1] = int(nodes.nodelist.index(templine[k2].split('=')[1]))
            if templine[k2].split('=')[0] == 'phases':
                lines.phases[k1] = templine[k2].split('=')[1]
                lines.PH[:,k1] = PHmat[:,phlist.index(templine[k2].split('=')[1])]
            if templine[k2].split('=')[0] == 'config':
                lines.config[k1] = templine[k2].split('=')[1]
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

    print(lines.TXnode)
    print(lines.RXnode)

    print(lines.TXnum)
    print(lines.RXnum)

    print(lines.phases)
    print(lines.PH)

    print(lines.config)
    print(lines.length)


    # Find the lines into and out of nodes

    intemp = np.zeros(nnode)
    outtemp = np.zeros(nnode)
    for k1 in range(0,nline):
        intemp[lines.RXnum[k1]]+=1
        outtemp[lines.TXnum[k1]]+=1
    print(intemp, int(np.amax(intemp)))
    print(outtemp, int(np.amax(outtemp)))

    nodes.inmat = -np.ones((int(np.amax(intemp)),nnode), dtype='int')
    nodes.outmat = -np.ones((int(np.amax(outtemp)),nnode), dtype='int')
    for k1 in range(0,nodes.nnode):
        incount = 0
        outcount = 0
        for k2 in range(0,lines.nline):
            if k1 == lines.RXnum[k2]:
                nodes.inmat[incount,k1] = k2
                incount+=1
            if k1 == lines.TXnum[k2]:
                nodes.outmat[outcount,k1] = k2
                outcount+=1

    print(nodes.inmat)
    print(nodes.outmat)


    # Line Config

    class configclass():
        pass
    configs = configclass()

    nconf = len(configlinearray)
    configs.nconf = nconf
    configs.conflist = [None]*nconf
    configs.FZpl = np.zeros((3,3*nline), dtype='complex')
    #print(a)

    for k1 in range(0,len(configlinearray)):
        templine = configlinearray[k1]

        tempFZpl = np.zeros((3,3))
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

        tempFZpl = np.array([[raa + 1j*xaa, rab + 1j*xab, rac + 1j*xac], \
                            [rab + 1j*xab, rbb + 1j*xbb, rbc + 1j*xbc], \
                            [rac + 1j*xac, rbc + 1j*xbc, rcc + 1j*xcc]])

        #print(k1, configs.conflist[k1], tempFZpl)

        for k2 in range(0,len(templine)):
            if templine[k2].split('=')[0] == 'unit':
                if templine[k2].split('=')[1] == 'ohm/m':
                    pass
                if templine[k2].split('=')[1] == 'ohm/km':
                    tempFZpl = tempFZpl/1000
                if templine[k2].split('=')[1] == 'ohm/ft' or templine[k2].split('=')[1] == 'ohm/feet':
                    tempFZpl = tempFZpl/0.3048
                if templine[k2].split('=')[1] == 'ohm/mi' or templine[k2].split('=')[1] == 'ohm/mile':
                    tempFZpl = tempFZpl/1609.34
        #print(tempFZpl)
        configs.FZpl[:,3*k1:3*(k1+1)] = tempFZpl

    #print(configs.FZpl)


    # Reconcile Lines and Configurations

    lines.FZ = np.zeros((3,3*nline), dtype='complex')
    lines.FZpu = np.zeros((3,3*nline), dtype='complex')
    lines.FY = np.zeros((3,3*nline), dtype='complex')
    lines.FYpu = np.zeros((3,3*nline), dtype='complex')

    for k1 in range(0,nline):
        confnum = configs.conflist.index(lines.config[k1])
        tempFZ = lines.length[k1]*configs.FZpl[:,3*confnum:3*(confnum+1)].reshape((3,3))

        for kph in range(0,3):
            if lines.PH[kph,k1] == 0:
                tempFZ[:,kph] = 0
                tempFZ[kph,:] = 0

        lines.FZ[:,3*k1:3*(k1+1)] = tempFZ
        tempFZpu = tempFZ/network.Zbase
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
        lines.FY = lines.FYpu/network.Zbase

    lines.FRpu = lines.FZpu.real
    lines.FXpu = lines.FZpu.imag

    lines.FGpu = lines.FYpu.real
    lines.FBpu = lines.FYpu.imag

    print('\nLINES')
    for k1 in range(0,nline):
        print(lines.FZpu[:,3*k1:3*(k1+1)])
    for k1 in range(0,nline):
        print(lines.FYpu[:,3*k1:3*(k1+1)])


    '''
    for k1 in range(0,nline):
        tempFZpu = lines.FZpu[:,3*k1:3*(k1+1)]
        tempFYpu = lines.FYpu[:,3*k1:3*(k1+1)]
        tempreal = (tempFZpu@tempFYpu).real
        tempimag = (tempFZpu@tempFYpu).imag
        tempreal[np.abs(tempreal) <= 1e-9] = 0
        tempimag[np.abs(tempimag) <= 1e-9] = 0
        print(tempreal + 1j*tempimag)
        #print(tempimag)
    '''


    # Loads

    class loadsclass():
        pass

    loads = loadsclass()

    loads.conn = []
    loads.p = np.zeros((3,nodes.nnode))
    loads.q = np.zeros((3,nodes.nnode))
    loads.aPQ = np.zeros((3,nodes.nnode))
    loads.aI = np.zeros((3,nodes.nnode))
    loads.aZ = np.zeros((3,nodes.nnode))

    for k1 in range(0,len(loadlinearray)):
        templine = loadlinearray[k1]

        for k2 in range(0,len(templine)):
            if templine[k2].split('=')[0] == 'nodename':
                knode = nodes.nodelist.index(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'nodenum':
                knode = int(templine[k2].split('=')[1])

            if templine[k2].split('=')[0] == 'conn':
                pass

            if templine[k2].split('=')[0] == 'phases':
                if templine[k2].split('=')[1] == 'a':
                    kph = 0
                elif templine[k2].split('=')[1] == 'b':
                    kph = 1
                elif templine[k2].split('=')[1] == 'c':
                    kph = 2

            if templine[k2].split('=')[0] == 'type':
                pass

            if templine[k2].split('=')[0] == 'apq':
                loads.aPQ[kph,knode] = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'ai':
                loads.aI[kph,knode] = float(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'az':
                loads.aZ[kph,knode] = float(templine[k2].split('=')[1])

            if templine[k2].split('=')[0] == 'real':
                loads.p[kph,knode] = float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[1] == 'W':
                    pass
                elif templine[k2+1].split('=')[1] == 'kW':
                    loads.p[kph,knode] = 1e3*loads.p[kph,knode]
                if templine[k2+1].split('=')[1] == 'MW':
                    loads.p[kph,knode] = 1e6*loads.p[kph,knode]

            if templine[k2].split('=')[0] == 'reac':
                loads.q[kph,knode] = float(templine[k2].split('=')[1])
                if templine[k2+1].split('=')[1] == 'VAr':
                    pass
                elif templine[k2+1].split('=')[1] == 'kVAr':
                    loads.q[kph,knode] = 1e3*loads.q[kph,knode]
                elif templine[k2+1].split('=')[1] == 'MVAr':
                    loads.q[kph,knode] = 1e6*loads.q[kph,knode]

    loads.ppu = loads.p/network.Sbase
    loads.qpu = loads.q/network.Sbase

    loads.s = loads.p + 1j*loads.q
    loads.spu = loads.ppu + 1j*loads.qpu

    print('\nLOADS')
    print(1.125*loads.s)
    print(1.125*loads.spu)
    

    # Capacitors

    class capsclass():
        pass
    caps = capsclass()

    caps.cap = np.zeros((3,nodes.nnode))
    caps.cappu = np.zeros((3,nodes.nnode))

    for k1 in range(0,len(capacitorlinearray)):
        templine = capacitorlinearray[k1]

        for k2 in range(0,len(templine)):
            if templine[k2].split('=')[0] == 'nodename':
                knode = nodes.nodelist.index(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'nodenum':
                knode = int(templine[k2].split('=')[1])

            if templine[k2].split('=')[0] == 'conn':
                if templine[k2].split('=')[1] == 'wye':
                    pass
                if templine[k2].split('=')[1] == 'delta':
                    pass

            if templine[k2].split('=')[0] == 'phase':
                if templine[k2].split('=')[1] == 'a':
                    kph = 0
                elif templine[k2].split('=')[1] == 'b':
                    kph = 1
                elif templine[k2].split('=')[1] == 'c':
                    kph = 2

            if templine[k2].split('=')[0] == 'reac':
                caps.cap[kph,knode] = float(templine[k2].split('=')[1])

            if templine[k2].split('=')[0] == 'unit':
                if templine[k2].split('=')[1] == 'VAr':
                    pass
                if templine[k2].split('=')[1] == 'kVAr':
                    caps.cap[kph,knode] = 1e3*caps.cap[kph,knode]

    caps.cappu = caps.cap/network.Sbase
    
    print('\nCAPS')
    print(caps.cap)
    print(caps.cappu)


    # Controllers

    '''
    class Controllers():

        def __init__(self):

            self.wmax = np.zeros((3,nodes.nnode))
            self.wmaxpu = np.zeros((3,nodes.nnode))

        def calculate(self):
            pass # do something

    controllers = Controllers()
    controllers.calculate()
    '''
    
    class controllersclass():
        pass
    
    controllers = controllersclass()
    controllers.wmax = np.zeros((3,nnode))
    controllers.wmaxpu = np.zeros((3,nnode))

    for k1 in range(0,len(controllerlinearray)):
        templine = controllerlinearray[k1]

        for k2 in range(0,len(templine)):
            if templine[k2].split('=')[0] == 'nodename':
                knode = nodes.nodelist.index(templine[k2].split('=')[1])
            if templine[k2].split('=')[0] == 'nodenum':
                knode = int(templine[k2].split('=')[1])

            if templine[k2].split('=')[0] == 'conn':
                if templine[k2].split('=')[1] == 'wye':
                    pass
                if templine[k2].split('=')[1] == 'delta':
                    pass

            if templine[k2].split('=')[0] == 'phases':
                if templine[k2].split('=')[1] == 'a':
                    kph = 0
                elif templine[k2].split('=')[1] == 'b':
                    kph = 1
                elif templine[k2].split('=')[1] == 'c':
                    kph = 2

            if templine[k2].split('=')[0] == 'mag':
                controllers.wmax[kph,knode] = float(templine[k2].split('=')[1])

            if templine[k2].split('=')[0] == 'unit':
                if templine[k2].split('=')[1] == 'VAr':
                    pass
                if templine[k2].split('=')[1] == 'kVAr':
                    controllers.wmax[kph,knode] = 1e3*controllers.wmax[kph,knode]
                if templine[k2].split('=')[1] == 'MVAr':
                    controllers.wmax[kph,knode] = 1e6*controllers.wmax[kph,knode]

    controllers.wmaxpu = controllers.wmax/network.Sbase
    
    print('\nCONTROLLERS')
    print(0.25*controllers.wmax)
    print(0.25*controllers.wmaxpu)
    
    

    return network, nodes, lines, configs, loads, caps, controllers