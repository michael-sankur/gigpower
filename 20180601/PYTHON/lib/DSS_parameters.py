import opendssdirect as dss
import re
import numpy as np
import time
from lib.zip_parameters import *
from lib.helper import transformer_regulator_parameters, nominal_load_values, nominal_cap_arr, wpu_final_arr
def relevant_openDSS_parameters(t, vvc_objects):

    dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[0])

    #TRANSFORMER, VOLTAGE REGULATOR lines
    tf_no = len(dss.Transformers.AllNames()) - len(dss.RegControls.AllNames()) #number of transformers
    tf_bus, vr_bus, _, _, _, vr_no, _ = transformer_regulator_parameters()

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
        pattern = r"(\w+)."
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
        TXnode[lineoutidx] = dss.Circuit.AllBusNames()[vr_bus[1, line]]
        RXnode[lineoutidx] = dss.Circuit.AllBusNames()[vr_bus[0, line]]
        TXnum[lineoutidx] = vr_bus[1, line]
        RXnum[lineoutidx] = vr_bus[0, line]


    #spu, apq, ai, az
    spu = np.zeros((3,nnode))
    ppu = np.zeros((3,nnode))
    qpu = np.zeros((3,nnode))
    aPQ = np.zeros((3,nnode))
    aI = np.zeros((3,nnode))
    aZ = np.zeros((3,nnode))

    beta_S = get_S()
    beta_I = get_I()
    beta_Z = get_Z()
       
    for n in range(len(dss.Loads.AllNames())): #go through the loads
        dss.Loads.Name(dss.Loads.AllNames()[n]) #set the load
        load_phases = [0, 0, 0] #instantiate load phases as all non-existent
        load_data = dss.CktElement.BusNames()[0].split('.')[1:]
        knode = dss.Circuit.AllBusNames().index((dss.CktElement.BusNames()[0].split('.')[0])) #match busname to index
        for i in load_data:
            phase = int(i)
            load_phases[phase-1] = 1
            aPQ[phase - 1, knode] = beta_S
            aZ[phase - 1, knode] = beta_Z
            aI[phase - 1, knode] = beta_I
        if len(load_data) == 0:
            load_phases = [1, 1, 1]  
            for p in range(len(load_phases)):            
                aPQ[p, knode] = beta_S
                aZ[p, knode] = beta_Z
                aI[p, knode] = beta_I    
                
    #ppu, qpu = load_values(t)
    ppu, qpu = nominal_load_values(-1)
    spu = (ppu + 1j * qpu)

    #cappu, wpu, vvcpu
    cappu = nominal_cap_arr()
 
    wpu = np.zeros((3, nnode))
    wpu = wpu_final_arr(nnode, vvc_objects)

    vvcpu = np.zeros((3,nnode))

    return TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu
