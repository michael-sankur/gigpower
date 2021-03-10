import opendssdirect as dss
import re
import numpy as np
import time
def relevant_openDSS_parameters(fn, t):

    dss.run_command('Redirect ' + fn)

    tf_no = len(dss.Transformers.AllNames()) - len(dss.RegControls.AllNames()) #number of transformers
    vr_no = len(dss.RegControls.AllNames()) #number of voltage regulators

    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])

    #BASE values
    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0
    Ibase = Sbase/Vbase
    #Zbase = Vbase/Ibase

    #TRANSFORMER, VOLTAGE REGULATOR lines
    tf_bus = np.zeros((2, tf_no), dtype = int) #tf has in and out bus
    vr_bus = np.zeros((3, vr_no), dtype = int) #vr has in and out bus and phase
    tf_count = 0
    vr_count = 0
    whichone = -1
    for tf in range(len(dss.Transformers.AllNames())):
        dss.Transformers.Name(dss.Transformers.AllNames()[tf])
        for i in range(2):
            if dss.Transformers.AllNames()[tf] in dss.RegControls.AllNames():
                bus = dss.CktElement.BusNames()[i].split('.')
                vr_bus[i, vr_count] = int(dss.Circuit.AllBusNames().index(bus[0]))
                vr_bus[2, vr_count] = int(bus[1]) #phase
                whichone = 0
            else:
                tf_bus[i, tf_count] =  int(dss.Circuit.AllBusNames().index(dss.CktElement.BusNames()[i])) #stuff the in and out bus of the tf into an array
                whichone = 1
        if whichone == 1:
            tf_count += 1
        else:
            vr_count += 1
        whichone = -1

    # vr_no = 1
    # vr_bus = vr_bus[0:2, 0:1]

    nline = len(dss.Lines.AllNames()) + tf_no + (2* vr_no)
    nnode = len(dss.Circuit.AllBusNames())

    line_in_idx_vr = range(len(dss.Lines.AllNames()) + tf_no, nline, 2)
    line_out_idx_vr = range(len(dss.Lines.AllNames()) + tf_no + 1, nline, 2)

    line_idx_tf = range(len(dss.Lines.AllNames()), len(dss.Lines.AllNames()) + tf_no)

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
        bus1 = dss.Lines.Bus1().split('.')[0]
        bus2 = dss.Lines.Bus2().split('.')[0]
        TXnode[line] = bus1
        RXnode[line] = bus2
        TXnum[line] = dss.Circuit.AllBusNames().index(TXnode[line])
        RXnum[line] = dss.Circuit.AllBusNames().index(RXnode[line])

    #TF
    for line in range(len(line_idx_tf)):
        print('is there a trnasformer?')
        lineidx = line_idx_tf[line]
        TXnode[lineidx] = dss.Circuit.AllBusNames()[tf_bus[0, line]] #bus name
        RXnode[lineidx] = dss.Circuit.AllBusNames()[tf_bus[1, line]]
        TXnum[lineidx] = tf_bus[0, line] #bus index
        RXnum[lineidx] = tf_bus[1, line]
    #VR in
    for line in range(len(line_in_idx_vr)):
        print('is there a voltage regulator?')
        lineinidx = line_in_idx_vr[line]

        TXnode[lineinidx] = dss.Circuit.AllBusNames()[vr_bus[0, line]]
        RXnode[lineinidx] = dss.Circuit.AllBusNames()[vr_bus[1, line]]
        TXnum[lineinidx] = vr_bus[0, line]
        RXnum[lineinidx] = vr_bus[1, line]
    #VR out
    for line in range(len(line_out_idx_vr)):
        print('is there a voltage regulator part 2?')
        lineoutidx = line_out_idx_vr[line]
        TXnode[lineoutidx] = dss.Circuit.AllBusNames()[vr_bus[0, line]]
        RXnode[lineoutidx] = dss.Circuit.AllBusNames()[vr_bus[1, line]]
        TXnum[lineoutidx] = vr_bus[0, line]
        RXnum[lineoutidx] = vr_bus[1, line]

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

    cap_dict()
    wpu = np.zeros((3, nnode))
    vvcpu = np.zeros((3,nnode))

    return TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu
