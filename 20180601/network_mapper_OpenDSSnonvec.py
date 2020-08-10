import numpy as np
import math
import opendssdirect as dss
import re

printflag = 1

def network_mapper_function(fn):



    baselinearray = [];
    nodelinearray = [];
    linelinearray = [];
    configlinearray = [];
    loadlinearray = [];
    capacitorlinearray = [];
    controllerlinearray = [];
    vvclinearray = [];

    dss.run_command('Redirect ' + str(fn))
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
    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])

    Vbase = dss.Bus.kVBase()
    print(dss.Bus.kVBase())
    Sbase = 1000000.0

    # Ibase = Sbase/Vbase
    # Zbase = Vbase/Ibase

    base.Vbase = Vbase
    base.Ibase = 360.8440864871107
    base.Zbase = 7.679992838399998
    base.Sbase = 1e6
    Ibase = base.Ibase
    Zbase = base.Zbase



    print('BASE')
    print('Vbase:', base.Vbase)
    print('Sbase:', base.Sbase)
    print('Ibase:', base.Ibase)
    print('Zbase:', base.Zbase)
    #########################
    # Nodes
    #########################

    # Number of nodes
    nnode = len(dss.Circuit.AllBusNames())
    nodes.nnode = nnode
    nodes.nodelist = [None]*nnode
    nodes.phases = [None]*nnode
    nodes.PH = np.zeros((3,nnode), dtype='int')

    # Iterate through node in configuration file. Nodes are added in order of their
    # their placement in the configuration file. The slacknode should be the
    # first node in the configuration file (for the time being)

    count = 0
    dictionary = {}
    for k2 in range(len(dss.Circuit.AllNodeNames())):
        for i in range(1, 4):
            pattern = r"\.%s" % (str(i))

            m = re.findall(pattern, dss.Circuit.AllNodeNames()[k2])
            a, b = dss.Circuit.AllNodeNames()[k2].split('.')
            if m and a in dictionary:
                temp = dictionary[a]
                temp[i - 1] = 1
                dictionary[a] = temp
            elif m and a not in dictionary:
                dictionary[a] = [0, 0, 0]
                temp = dictionary[a]
                temp[i - 1] = 1
                dictionary[a] = temp
    for key, value in dictionary.items():
        nodes.nodelist[count] = key
        phase_str = ""
        if value[0] == 1:
            phase_str = phase_str + "a"
        if value[1] == 1:
            phase_str = phase_str + "b"
        if value[2] == 1:
            phase_str = phase_str + "c"

        nodes.PH[:, count] = value
        nodes.phases[count] = phase_str

        count = count + 1

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
    nline = len(dss.Lines.AllNames())
    lines.nline = nline
    lines.TXnode = [None]*nline #name of incoming line's bus
    lines.RXnode = [None]*nline #name of outgoing bus on line
    lines.TXnum = np.zeros((nline), dtype='int') #int value, do as dict
    lines.RXnum = np.zeros((nline), dtype='int') #int value
    lines.phases = [None]*nline #phase of the line
    lines.PH = np.zeros((3,nline), dtype='int') #ph of the line
    lines.config = [None]*nline
    lines.length = np.zeros((nline))

    busIdx = {}
    for i in range(len(dss.Circuit.AllBusNames())):
        busIdx[dss.Circuit.AllBusNames()[i]] = i

    print(busIdx)
    for line in range(len(dss.Lines.AllNames())):
        dss.Lines.Name(dss.Lines.AllNames()[line]) #set the line
        bus1 = dss.Lines.Bus1()
        bus2 = dss.Lines.Bus2()
        pattern = r"(\w+)."

        lines.TXnode[line] = re.findall(pattern, bus1)[0]
        lines.RXnode[line] = re.findall(pattern, bus2)[0]
        lines.TXnum[line] = busIdx[lines.TXnode[line]]
        lines.RXnum[line] = busIdx[lines.RXnode[line]]
        lines.length[line] = dss.Lines.Length()
        lines.config[line] = dss.Lines.LineCode()

    line_phase_dict = {}
    for k2 in range(len(dss.Lines.AllNames())):
        for i in range(1, 4):
            pattern = r"\.%s" % (str(i))
            dss.Lines.Name(dss.Lines.AllNames()[k2])
            m = re.findall(pattern, dss.Lines.Bus1())
            if m and dss.Lines.AllNames()[k2] in line_phase_dict :
                temp = line_phase_dict [dss.Lines.AllNames()[k2]]
                temp[i - 1] = 1
                line_phase_dict [dss.Lines.AllNames()[k2]] = temp
            elif m and dss.Lines.AllNames()[k2] not in line_phase_dict :
                line_phase_dict [dss.Lines.AllNames()[k2]] = [0, 0, 0]
                temp = line_phase_dict [dss.Lines.AllNames()[k2]]
                temp[i - 1] = 1
                line_phase_dict [dss.Lines.AllNames()[k2]] = temp
    count = 0
    for key, value in line_phase_dict.items():
        lines.PH[:, count] = value
        phase_str = ""
        if value[0] == 1:
            phase_str = phase_str + "a"
        if value[1] == 1:
            phase_str = phase_str + "b"
        if value[2] == 1:
            phase_str = phase_str + "c"
        print(value)

        lines.phases[count] = phase_str
        count += 1

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


    configs.nconf = len(lines.config)
    nconf = configs.nconf
    configs.conflist = [None]*nconf
    # Impedance matrices for configs [pu/m]
    configs.FZpupl = np.zeros((3,3*nline), dtype='complex')

    print(lines.PH)
    # Iterate through line configurations
    for k1 in range(0,len(lines.config)):
        dss.Lines.Name(dss.Lines.AllNames()[k1])
        linecode = dss.Lines.LineCode()
        rmat = dss.Lines.RMatrix()
        xmat = dss.Lines.XMatrix()

        rtemp = np.zeros((3, 3))
        xtemp = np.zeros((3, 3))

        if np.all(lines.PH[:, k1] ==  [1, 1, 0]):
            rtemp[0, 0] = rmat[0]
            rtemp[0, 1] = rmat[1]
            rtemp[1, 0] = rmat[2]
            rtemp[1, 1] = rmat[3]
            xtemp[0, 0] = xmat[0]
            xtemp[0, 1] = xmat[1]
            xtemp[1, 0] = xmat[2]
            xtemp[1, 1] = xmat[3]
        elif np.all(lines.PH[:, k1]==[1, 0, 1]):
            rtemp[0, 0] = rmat[0]
            rtemp[2, 0 ] = rmat[2]
            rtemp[0, 2] = rmat[1]
            rtemp[2, 2] = rmat[3]
            xtemp[0, 0] = xmat[0]
            xtemp[2, 0] = xmat[2]
            xtemp[0, 2] = xmat[1]
            xtemp[2, 2] = xmat[3]
        elif np.all(lines.PH[:, k1]== [0,1, 1]):
            rtemp[1, 1] = rmat[0]
            rtemp[1, 2 ] = rmat[2]
            rtemp[2, 1] = rmat[1]
            rtemp[2, 2] = rmat[3]
            xtemp[1, 1] = xmat[0]
            xtemp[1, 2 ] = xmat[2]
            xtemp[2, 1] = xmat[1]
            xtemp[2, 2] = xmat[3]
        elif np.all(lines.PH[:,k1]== [1, 1, 1]):
            rtemp = np.reshape(rmat, (3, 3))
            xtemp = np.reshape(xmat,(3, 3))
        elif np.all(lines.PH[:,k1]== [1, 0, 0]):
            rtemp[0, 0] = rmat[0]
            xtemp[0, 0] = xmat[0]
        elif np.all(lines.PH[:,k1] == [0, 1, 0]):
            rtemp[1, 1] = rmat[0]
            xtemp[1, 1] = xmat[0]
        else:
            rtemp[2, 2] = rmat[0]
            xtemp[2, 2] = xmat[0]
        # Impedance (resistance and reactance) values for current line configuration
        tempFZpupl = np.zeros((3,3))
        tempFZpupl = rtemp + 1j*xtemp

        # 3x3 impedance per unit length matrix for current line config [pu/m]
        configs.FZpupl[:,3*k1:3*(k1+1)] = tempFZpupl/Zbase
        print(Zbase)
    configs.conflist = lines.config

    if printflag == 1:
        print()
        print('CONFIGS')
        print('nconf:', configs.nconf)
        print('conflist:', configs.conflist)
        for k1 in range(configs.nconf):
            print(configs.conflist[k1], '- FZpl:\n', configs.FZpupl[:,3*k1:3*(k1+1)])

    # Reconcile Lines and Configurations

    # Impednace matrices for line [pu]
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
    for k in range(len(dss.Loads.AllNames())):
        dss.Loads.Name(dss.Loads.AllNames()[k])
        no_pattern = r"_([\w]+?)_"
        knode = re.findall(no_pattern, dss.Loads.Name())
        knode = nodes.nodelist.index(knode[0])
        ph_pattern = r"_([\w]?)_"
        kph = re.findall(ph_pattern, dss.Loads.Name())

        if kph[0] == 'a':
            kph = 0
        elif kph[0] == 'b':
            kph = 1
        elif kph[0] == 'c':
            kph = 2
        loads.ppu[kph,knode] = dss.Loads.kW()
        loads.qpu[kph,knode] = dss.Loads.kvar()

    apq = 0.85
    ai = 0.05
    az = 0.1

        # if templine[k2].split('=')[0] == 'nodenum':
        #     knode = int(templine[k2].split('=')[1])

    # Complex loads [pu]
    loads.spu = loads.ppu + 1j*loads.qpu
    loads.spu_nominal = loads.ppu +1j*loads.qpu

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

    # # Iterate through capacitors in configuration file
    # for k1 in range(0,len(dss.Capacitors.AllNames())):
    #     dss.Capacitors.Name(dss.Capacitors.AllNames()[k1])
    #     knode = dss.Capacitors.Name()
    #         # if templine[k2].split('=')[0] == 'nodenum':
    #         #     knode = int(templine[k2].split('=')[1])
    #         # # Capacitor phase(s)
    #         # if templine[k2].split('=')[0] == 'phase':
    #         #     if templine[k2].split('=')[1] == 'a':
    #         #         kph = 0
    #         #     elif templine[k2].split('=')[1] == 'b':
    #         #         kph = 1
    #         #     elif templine[k2].split('=')[1] == 'c':
    #         #         kph = 2
    #
    #     caps.cappu[kph, knode] = dss.Capacitors.kvar()

    if printflag == 1:
        print()
        print('CAPS')
        print('cappu:\n', caps.cappu)

    #########################
    # Controllers
    #########################

    # Controller dispatch [pu]
    cons.wpu = np.zeros((3,nodes.nnode))

    if printflag == 1:
        print()
        print('CONS')
        print('wpu:\n', cons.wpu)

    # VVC dispatch [pu]
    vvc.vvcpu = np.zeros((3,nnode))


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
