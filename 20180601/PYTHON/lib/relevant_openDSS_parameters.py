import opendssdirect as dss
import re
import numpy as np
def relevant_openDSS_parameters(fn):

    dss.run_command('Redirect ' + fn)
    dss.Solution.Solve()
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())
    nodelist = [None]*nnode

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


    for k in range(len(dss.Loads.AllNames())):
        dss.Loads.Name(dss.Loads.AllNames()[k])
        no_pattern = r"_([\w]+?)_"
        knode = re.findall(no_pattern, dss.Loads.Name())
        knode = nodelist.index(knode[0])
        ph_pattern = r"_([\w]?)_"
        kph = re.findall(ph_pattern, dss.Loads.Name())

        if kph[0] == 'a':
            kph = 0
        elif kph[0] == 'b':
            kph = 1
        elif kph[0] == 'c':
            kph = 2

        aPQ[kph, knode] = 1 #temporary
        aI[kph,knode] = 0.00
        aZ[kph, knode] = 0.00
        ppu[kph,knode] = dss.Loads.kW() * 1000 / 1000000.0
        qpu[kph,knode] = dss.Loads.kvar() * 1000 / 1000000.0
    spu = ppu + 1j * qpu

    #cappu, wpu, vvcpu
    cappu = np.zeros((3,nnode))
    wpu = np.zeros((3, nnode))
    vvcpu = np.zeros((3,nnode))

    return TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu
