import numpy as np
import math
import opendssdirect as dss
import re

printflag = 1

def network_mapper_function(fn, t):
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
    #print(dss.Circuit.AllBusNames())
    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])

    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0

    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase

    base.Vbase = Vbase
    base.Ibase = Ibase
    base.Zbase = Zbase
    base.Sbase = Sbase

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

    for line in range(len(dss.Lines.AllNames())):
        dss.Lines.Name(dss.Lines.AllNames()[line]) #set the line
        bus1 = dss.Lines.Bus1()
        bus2 = dss.Lines.Bus2()
        pattern = r"(\w+)."

        lines.TXnode[line] = re.findall(pattern, bus1)[0]
        lines.RXnode[line] = re.findall(pattern, bus2)[0]
        lines.TXnum[line] = busIdx[lines.TXnode[line]]
        lines.RXnum[line] = busIdx[lines.RXnode[line]]
        lines.length[line] = dss.Lines.Length() #* 0.3048
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
                line_phase_dict[dss.Lines.AllNames()[k2]] = [0, 0, 0]
                temp = line_phase_dict[dss.Lines.AllNames()[k2]]
                temp[i - 1] = 1
                line_phase_dict[dss.Lines.AllNames()[k2]] = temp
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

    # Iterate through line configurations
    for k1 in range(0,len(lines.config)):
        dss.Lines.Name(dss.Lines.AllNames()[k1])
        linecode = dss.Lines.LineCode()
        dss.LineCodes.Name(linecode)
        rmat = dss.LineCodes.Rmatrix()
        xmat = dss.LineCodes.Xmatrix()

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
            rtemp[2, 0] = rmat[2]
            rtemp[0, 2] = rmat[1]
            rtemp[2, 2] = rmat[3]
            xtemp[0, 0] = xmat[0]
            xtemp[2, 0] = xmat[2]
            xtemp[0, 2] = xmat[1]
            xtemp[2, 2] = xmat[3]
        elif np.all(lines.PH[:, k1]== [0,1, 1]):
            rtemp[1, 1] = rmat[0]
            rtemp[1, 2] = rmat[2]
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
        configs.FZpupl[:,3*k1:3*(k1+1)] = tempFZpupl/Zbase#/1609.34
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
        dss.Lines.Name(dss.Lines.AllNames()[k1])
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

    def get_bus_idx(bus):
        pattern =  r"(\w+)\."
        bus_ph_free = re.findall(pattern, bus)
        k =  dss.Circuit.AllBusNames().index(bus_ph_free[0])
        return k

    # Iterate through loads in configuration
    # accomodate multiphase and multiple loads at 1 bus
    def load_order_f(): #dict {bus: no_loads}
        load_order = {}
        for n in range(len(dss.Loads.AllNames())):
            dss.Loads.Name(dss.Loads.AllNames()[n])
            pattern =  r"(\w+)\."
            load_bus = re.findall(pattern, dss.CktElement.BusNames()[0])
            if load_bus[0] not in load_order: #if bus DNE in dict
                load_order[load_bus[0]] = 1
            elif load_bus[0] in load_order: #if bus E in dict
                load_order[load_bus[0]] += 1 #add 1 onto count of loads
        return load_order
    load_order_list = load_order_f()

    def load_values():
        load_ph_arr = np.zeros((nnode, max(load_order_list.values()), 3))
        # load_kw_arr = np.zeros((nnode, max(load_order_list.values())))
        # load_kvar_arr = np.zeros((nnode, max(load_order_list.values())))
        load_kw_arr_ph = np.zeros((3, nnode))
        load_kvar_arr_ph = np.zeros((3, nnode))
        loads_apq = np.zeros((3, nnode))
        loads_az = np.zeros((3,nnode))

        for load in range(len(dss.Loads.AllNames())):
            dss.Loads.Name(dss.Loads.AllNames()[load])
            pattern =  r"(\w+)\."
            load_bus = re.findall(pattern, dss.CktElement.BusNames()[0]) #bus name for load
            load_ph_arr_temp = [0, 0, 0] #ph of load
            for i in range(1, 4): #update existing ph of load
                pattern = r"\.%s" % (str(i))
                load_ph = re.findall(pattern, dss.CktElement.BusNames()[0])
                if load_ph:
                    load_ph_arr_temp[i - 1] = 1
            for j in range(max(load_order_list.values())):
                idxbs = dss.Circuit.AllBusNames().index(load_bus[0]) #what index is bus at in allbusnames()
                if np.all(load_ph_arr[idxbs, j,:] == [0, 0, 0]): #if row is empty
                    load_ph_arr[idxbs, j, :] = load_ph_arr_temp #place ph array there
                    # load_kw_arr[idxbs, j] = dss.Loads.kW() / sum(load_ph_arr_temp)
                    # load_kvar_arr[idxbs, j] = dss.Loads.kvar() / sum(load_ph_arr_temp)
                    for i in range(len(load_ph_arr_temp)):
                        if load_ph_arr_temp[i] == 1:
                            load_kw_arr_ph[i,idxbs] += (dss.Loads.kW() * 1e3 / Sbase / sum(load_ph_arr_temp))
                            load_kvar_arr_ph[i,idxbs] += (dss.Loads.kvar() * 1e3 / Sbase / sum(load_ph_arr_temp))
                            loads_apq[i,idxbs] = 1.0
                            loads_az[i, idxbs] = 0.0
                    break
        #return load_ph_arr, load_kw_arr, load_kvar_arr, load_kw_arr_ph, load_kvar_arr_ph, loads_apq, loads_az
        return load_kw_arr_ph, load_kvar_arr_ph, loads_apq, loads_az
    #load_ph_arr, load_kw_arr, load_kvar_arr,
    load_kw_arr_ph, load_kvar_arr_ph,loads_apq, loads_az = load_values()

    loads.ppu = load_kw_arr_ph
    loads.qpu = load_kvar_arr_ph
    loads.aPQ = loads_apq
    loads.aZ = loads_az

    #Complex loads [pu]

    if t == -1:
        var = 1
    else:
        var = (1 + 0.1*np.sin(2*np.pi*0.01*t))
    loads.spu = (loads.ppu + 1j*loads.qpu)*var
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

    # def cap_dict():
    #     cap_dict_kvar = {}
    #     cap_dict_kv = {}
    #     cap_ph_dict = {}
    #     for n in range(len(dss.Capacitors.AllNames())):
    #         dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
    #         pattern =  r"(\w+)\."  #if it is, is the phase present?
    #         m = re.findall(pattern, dss.CktElement.BusNames()[0])
    #         cap_dict_kvar[m[0]] = dss.Capacitors.kvar()
    #         cap_dict_kv[m[0]] = dss.Capacitors.kV()
    #         load_phases = [0, 0, 0]
    #         print(m)
    #         for i in range(1, 4): #if the phase is present, what other phases are
    #             pattern = r"\.%s" % (str(i))
    #             m2 = re.findall(pattern, dss.CktElement.BusNames()[0])
    #             print(m2)
    #             if m2:
    #                  if m[0] not in cap_ph_dict.keys():
    #                     cap_ph_dict[m[0]] = [0, 0, 0]
    #                     cap_ph_dict[m[0]][i-1] = 1
    #                  else:
    #                     cap_ph_dict[m[0]][i-1] = 1
    #         #cap_ph_dict[m[0]] = load_phases
    #     return cap_dict_kvar, cap_dict_kv, cap_ph_dict
    # cap_dict_kvar, cap_dict_kv, cap_ph_dict = cap_dict()
    # print(cap_dict_kvar)
    # print(cap_ph_dict)
    caps.cappu = np.zeros((3,nodes.nnode))

    # # # Iterate through capacitors in configuration file
    # #print(dss.Circuit.AllBusNames())
    # for ph in range(0, 3):
    #     for bus in range(len(dss.Circuit.AllBusNames())):
    #         if dss.Circuit.AllBusNames()[bus] in cap_dict_kvar.keys():
    #             #print(bus)
    #             if cap_ph_dict[dss.Circuit.AllBusNames()[bus]][ph] == 1:
    #                 caps.cappu[ph, bus] = cap_dict_kvar[dss.Circuit.AllBusNames()[bus]] * 1000 / Sbase / sum(cap_ph_dict[dss.Circuit.AllBusNames()[bus]])
    # #print(caps.cappu)

#-------------
    cap_dict = {} #dict {bus_name: kvar}
    cap_ph_dict = {} #dict {bus_name:1x3 ph array}
    for capname in range(len(dss.Capacitors.AllNames())):
        dss.Capacitors.Name(dss.Capacitors.AllNames()[capname])
        pattern =  r"(\w+)\."
        load_phases = [0, 0, 0] #what phases exist in capacitor
        for i in range(1, 4): #if the phase is present, what other phases are
            pattern = r"\.%s" % (str(i))
            m2 = re.findall(pattern, dss.CktElement.BusNames()[0])
            if m2:
                load_phases[i - 1] = 1 #update phase array
        cap_ph_dict[dss.CktElement.BusNames()[0]] = load_phases
        cap_dict[dss.CktElement.BusNames()[0]] = dss.Capacitors.kvar() * 1e3 / Sbase / sum(load_phases)

    for bus in cap_ph_dict.keys(): #update caps.cappu
        pattern =  r"(\w+)\."
        cap_bus = re.findall(pattern, bus)
        idxbs = dss.Circuit.AllBusNames().index(cap_bus[0]) #locate bus idx in allbusnames()
        for ph in range(3):
            if cap_ph_dict[bus][ph] == 1: #if the ph at bus == 1
                caps.cappu[ph, idxbs] += cap_dict[bus]

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
