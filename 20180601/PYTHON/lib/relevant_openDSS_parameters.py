import opendssdirect as dss
import re
import numpy as np
import time
def relevant_openDSS_parameters(fn, t):

    dss.run_command('Redirect ' + fn)
    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])

    #BASE values
    Vbase = dss.Bus.kVBase() * 1000
    Sbase = 1000000.0
    #Ibase = Sbase/Vbase
    #Zbase = Vbase/Ibase

    #TRANSFORMER, VOLTAGE REGULATOR lines
    tf_no = len(dss.Transformers.AllNames()) - len(dss.RegControls.AllNames()) #number of transformers
    vr_no = len(dss.RegControls.AllNames()) #number of voltage regulators
  
    tf_bus = np.zeros((5, tf_no), dtype = int) #r1 = in bus, r2 = out bus, r3-5 phases; c = tf
    vr_bus = np.zeros((5, vr_no), dtype = int) 
   
    vr_count = 0 
    vr_lines = 0
    tf_count = 0
    tf_lines = 0

    for vr in range(len(dss.RegControls.AllNames())):
        dss.RegControls.Name(dss.RegControls.AllNames()[vr])  
        for i in range(2): #start / end bus
            dss.Transformers.Name(dss.RegControls.Transformer())     
            bus = dss.CktElement.BusNames()[i].split('.')
            vr_bus[i, vr_count] = int(dss.Circuit.AllBusNames().index(bus[0])) 
        for n in range(len(bus[1:])):                 
            vr_lines += 1                    
            vr_bus[int(bus[1:][n]) + 1, vr_count] = int(bus[1:][n]) # phases that exist                   
        if len(bus) == 1: # unspecified phases, assume all 3 exist
            for n in range(1,4): 
                vr_lines += 1
                vr_bus[n+1, vr_count] = n   
        vr_count += 1 
   
    tf_bus_temp = np.zeros((2, 1))
    for tf in range(len(dss.Transformers.AllNames())):
        dss.Transformers.Name(dss.Transformers.AllNames()[tf]) 
        for i in range(2):     
            bus = dss.CktElement.BusNames()[i].split('.')        
            tf_bus_temp[i] = int(dss.Circuit.AllBusNames().index(bus[0])) 
            # stuff the in and out bus of the tf into an array  
        if not np.size(np.where(vr_bus[0, :] == tf_bus_temp[0])) == 0 and \
        not np.size(np.where(vr_bus[1, :] == tf_bus_temp[1])) == 0:     
            continue #if you have already seen the transformer from regulators, skip       
        
        tf_bus[0:2, tf_count] = tf_bus_temp[:, 0] #otherwise add to the tf_bus matrix
        for n in range(len(bus[1:])):                 
            tf_lines += 1                             
            tf_bus[int(bus[1:][n]) + 1, tf_count] = int(bus[1:][n])   
        if len(bus) == 1:
            for k in range(1,4): #if all phases are assumed to exist
                tf_lines += 1
                tf_bus[k+1, tf_count] = k             
        tf_count += 1

    nline = len(dss.Lines.AllNames()) + tf_no + (2*vr_no)
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
        bus1 = dss.Lines.Bus1()
        bus2 = dss.Lines.Bus2()
        pattern = r"(\w+)." #this appears to wrok
        try:
            TXnode[line] = re.findall(pattern, bus1)[0]
            RXnode[line] = re.findall(pattern, bus2)[0]
            TXnum[line] = dss.Circuit.AllBusNames().index(TXnode[line])
            RXnum[line] = dss.Circuit.AllBusNames().index(RXnode[line])
        except:
            pattern = r"(\w+)"
            TXnode[line] = re.findall(pattern, bus1)[0]
            RXnode[line] = re.findall(pattern, bus2)[0]
            TXnum[line] = dss.Circuit.AllBusNames().index(TXnode[line])
            RXnum[line] = dss.Circuit.AllBusNames().index(RXnode[line])

    #TF
    for line in range(len(line_idx_tf)):
        lineidx = line_idx_tf[line]
        TXnode[lineidx] = dss.Circuit.AllBusNames()[tf_bus[0, line]] #bus name
        RXnode[lineidx] = dss.Circuit.AllBusNames()[tf_bus[1, line]]
        TXnum[lineidx] = tf_bus[0, line] #bus index
        RXnum[lineidx] = tf_bus[1, line]
    #VR in
    for line in range(len(line_in_idx_vr)):
        lineinidx = line_in_idx_vr[line]
        TXnode[lineinidx] = dss.Circuit.AllBusNames()[vr_bus[0, line]]
        RXnode[lineinidx] = dss.Circuit.AllBusNames()[vr_bus[1, line]]
        TXnum[lineinidx] = vr_bus[0, line]
        RXnum[lineinidx] = vr_bus[1, line]
    #VR out
    for line in range(len(line_out_idx_vr)):
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
        load_data = dss.CktElement.BusNames()[0].split('.')[1:]
        knode = dss.Circuit.AllBusNames().index((dss.CktElement.BusNames()[0].split('.')[0])) #match busname to index
        for i in load_data:
            phase = int(i)
            load_phases[phase-1] = 1
        if len(load_data) == 0:
            load_phases = [1, 1, 1]        
        realstuff = dss.CktElement.Powers()[::2]
        imagstuff = dss.CktElement.Powers()[1::2]     
        rs = 0
        for ph in range(len(load_phases)):      
            if load_phases[ph] == 1:
                aPQ[ph, knode] = 1
                aZ[ph, knode] = 0
                # print(ph,dss.Circuit.AllBusNames()[knode])
                # print(realstuff[rs])
                # print("\n")
                ppu[ph, knode] += realstuff[rs] * 1e3 * var / Sbase 
                qpu[ph, knode] += imagstuff[rs] * 1e3 * var / Sbase
                rs += 1

    spu = (ppu + 1j * qpu)

    #cappu, wpu, vvcpu
    cappu = np.zeros((3,nnode))
 
    for n in range(len(dss.Capacitors.AllNames())):
        dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
        cap_data = dss.CktElement.BusNames()[0].split('.')

        idxbs = dss.Circuit.AllBusNames().index(cap_data[0])
        for ph in range(1, len(cap_data)):        
            cappu[int(cap_data[ph]) - 1, idxbs] += dss.Capacitors.kvar() * 1e3 / Sbase / (len(cap_data) - 1)


    wpu = np.zeros((3, nnode))
    vvcpu = np.zeros((3,nnode))

    return TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu
