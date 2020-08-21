import opendssdirect as dss
import re
import numpy as np
def relevant_openDSS_parameters(fn):

    dss.run_command('Redirect ' + fn)
    #dss.Solution.Solve()
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())

    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])

    #BASE values
    Vbase = dss.Bus.kVBase() * 1000  #@mike edit
    Sbase = 1000000.0  #@mike edit
    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase

    nodelist = [None]*nnode
    #NODE indexing
    TXnum = np.zeros((nline), dtype='int') #int value, do as dict
    RXnum = np.zeros((nline), dtype='int') #int value
    TXnode = [None]*nline #name of incoming line's bus
    RXnode = [None]*nline #name of outgoing bus on line
    #PH = np.zeros((3,nline), dtype='int')
    PH = np.zeros((3,nnode), dtype='int')
    #bus indexes
    busIdx = {}
    for i in range(len(dss.Circuit.AllBusNames())):
        busIdx[dss.Circuit.AllBusNames()[i]] = i

    #PH
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
    count = 0

    for key, value in dictionary.items():
        nodelist[count] = key
        PH[:, count] = value
        count += 1

    #TXnum and RXnum
    for line in range(len(dss.Lines.AllNames())):
        dss.Lines.Name(dss.Lines.AllNames()[line]) #set the line
        bus1 = dss.Lines.Bus1()
        bus2 = dss.Lines.Bus2()
        pattern = r"(\w+)."

        TXnode[line] = re.findall(pattern, bus1)[0]
        RXnode[line] = re.findall(pattern, bus2)[0]
        TXnum[line] = busIdx[TXnode[line]]
        RXnum[line] = busIdx[RXnode[line]]
    #spu, apq, ai, az
    spu = np.zeros((3,nnode))
    ppu = np.zeros((3,nnode))
    qpu = np.zeros((3,nnode))
    aPQ = np.zeros((3,nnode))
    aI = np.zeros((3,nnode))
    aZ = np.zeros((3,nnode))

    def get_bus_idx(bus):
        k = -1
        for n in range(len(dss.Circuit.AllBusNames())): #iterates over all the buses to see which index corresponds to bus
            if dss.Circuit.AllBusNames()[n] in bus:
                k = n
        return k

    for kph in range(0, 3):
        for k in range(len(dss.Circuit.AllBusNames())):
            for n in range(len(dss.Loads.AllNames())): #go through the loads
                dss.Loads.Name(dss.Loads.AllNames()[n]) #set the load
                if dss.Circuit.AllBusNames()[k] in dss.CktElement.BusNames()[0]: #check is the busname in the busname of the load
                    pattern =  r"\.%s" % (str(kph + 1)) #if it is, is the phase present?
                    m = re.findall(pattern, dss.CktElement.BusNames()[0])
                    if m:
                        load_phases = [0, 0, 0]
                        for i in range(1, 4): #if the phase is present, what other phases are
                            pattern = r"\.%s" % (str(i))
                            m2 = re.findall(pattern, dss.CktElement.BusNames()[0])
                            if m2:
                                load_phases[i - 1] = 1
                        knode = get_bus_idx(dss.Circuit.AllBusNames()[k])
                        if sum(load_phases) == 1:

                            aPQ[kph, knode] = 0.85 #temporary
                            aZ[kph,knode] = .15
                            ppu[kph,knode] = dss.Loads.kW() * 1e3 / Sbase
                            qpu[kph,knode] = dss.Loads.kvar() * 1e3 / Sbase
                        else:
                            aPQ[kph, knode] = 0.85 #temporary
                            aZ[kph,knode] = .15
                            ppu[kph,knode] =( dss.Loads.kW() + dss.Loads.kvar())* 1e3 / Sbase / sum(load_phases)
                            qpu[kph,knode] = ( dss.Loads.kW() + dss.Loads.kvar())* 1e3 / Sbase / sum(load_phases)

    spu = ppu + 1j * qpu

    #cappu, wpu, vvcpu
    cappu = np.zeros((3,nnode))
    wpu = np.zeros((3, nnode))
    vvcpu = np.zeros((3,nnode))

    return TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu
