import opendssdirect as dss
import re
import numpy as np
import time
def relevant_openDSS_parameters(fn, t):

    dss.run_command('Redirect ' + fn)
    #dss.Solution.Solve()
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())

    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])

    #BASE values
    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0
    Ibase = Sbase/Vbase
    Zbase = Vbase/Ibase

    nodelist = dss.Circuit.AllBusNames()
    #NODE indexing
    TXnum = np.zeros((nline), dtype='int')
    RXnum = np.zeros((nline), dtype='int')
    TXnode = [None]*nline
    RXnode = [None]*nline #name of outgoing bus on line
    PH = np.zeros((3,nnode), dtype='int')

    #PH
    for k2 in range(len(dss.Circuit.AllNodeNames())):
        a, b = dss.Circuit.AllNodeNames()[k2].split('.')
        PH[int(b) - 1, dss.Circuit.AllBusNames().index(a)] = 1

    #TXnum and RXnum
    for line in range(len(dss.Lines.AllNames())):
        dss.Lines.Name(dss.Lines.AllNames()[line]) #set the line
        bus1 = dss.Lines.Bus1()
        bus2 = dss.Lines.Bus2()
        pattern = r"(\w+)." #this appears to wrok

        TXnode[line] = re.findall(pattern, bus1)[0]
        RXnode[line] = re.findall(pattern, bus2)[0]
        TXnum[line] = dss.Circuit.AllBusNames().index(TXnode[line])
        RXnum[line] = dss.Circuit.AllBusNames().index(RXnode[line])

    #spu, apq, ai, az
    spu = np.zeros((3,nnode))
    ppu = np.zeros((3,nnode))
    qpu = np.zeros((3,nnode))
    aPQ = np.zeros((3,nnode))
    aI = np.zeros((3,nnode))
    aZ = np.zeros((3,nnode))

    if t == -1:
        var = 1
    else:
        var = (1 + 0.1*np.sin(2*np.pi*0.01*t))

    for n in range(len(dss.Loads.AllNames())): #go through the loads
        dss.Loads.Name(dss.Loads.AllNames()[n]) #set the load
        load_phases = [0, 0, 0] #instantiate load phases as all non-existent
        pattern =  r"(\w+)\."
        load_bus = re.findall(pattern, dss.CktElement.BusNames()[0]) #determine bus name
        knode = dss.Circuit.AllBusNames().index(load_bus[0]) #match busname to index
        for ph in range(0, 3):
            pattern =  r"\.%s" % (str(ph + 1))
            m = re.findall(pattern, dss.CktElement.BusNames()[0])
            if m: #if phase exists for load
                load_phases[ph - 1] = 1
                aPQ[ph, knode] = 1
                aZ[ph, knode] = 0
                ppu[ph, knode] = dss.Loads.kW()* 1e3 * var / Sbase #check these lines later
                qpu[ph, knode] = dss.Loads.kvar() * 1e3 * var / Sbase
        if sum(load_phases) > 1: #if it is a multiphase load
            for ph in range(0,3):
                ppu[ph, knode]/= sum(load_phases)
                qpu[ph, knode] /= sum(load_phases)
    spu = (ppu + 1j * qpu)

    #cappu, wpu, vvcpu
    cappu = np.zeros((3,nnode))

    def cap_dict():
        for n in range(len(dss.Capacitors.AllNames())):
            dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
            cap_data = dss.CktElement.BusNames()[0].split('.')

            idxbs = dss.Circuit.AllBusNames().index(cap_data[0])
            for ph in range(1, len(cap_data)):
                cappu[int(cap_data[ph]) - 1, idxbs] += dss.Capacitors.kvar() * 1e3 / Sbase / (len(cap_data) - 1)


    wpu = np.zeros((3, nnode))
    vvcpu = np.zeros((3,nnode))

    return TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu
