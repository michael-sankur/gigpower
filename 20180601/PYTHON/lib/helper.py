import numpy as np
import opendssdirect as dss
import re

def load_values(t):
    if t == -1:
        var = 1
    else:
        var = (1 + 0.1*np.sin(2*np.pi*0.01*t))
    
    nnode = len(dss.Circuit.AllBusNames())
    Sbase = 1000000.0

    ppu = np.zeros((3,nnode))
    qpu = np.zeros((3,nnode))

    for n in range(len(dss.Loads.AllNames())): #go through the loads
        dss.Loads.Name(dss.Loads.AllNames()[n]) #set the load
        load_phases = [0, 0, 0] #instantiate load phases as all non-existent
        load_data = dss.CktElement.BusNames()[0].split('.')[1:]
        knode = dss.Circuit.AllBusNames().index((dss.CktElement.BusNames()[0].split('.')[0])) #match busname to index
        for p in load_data:
            phase = int(p)
            load_phases[phase-1] = 1
        if len(load_data) == 0:
            load_phases = [1, 1, 1]        
        realstuff = dss.CktElement.Powers()[::2]
        imagstuff = dss.CktElement.Powers()[1::2]     
        rs = 0
        for ph in range(len(load_phases)):      
            if load_phases[ph] == 1:
                ppu[ph, knode] += realstuff[rs] * 1e3 * var / Sbase 
                qpu[ph, knode] += imagstuff[rs] * 1e3 * var / Sbase
                rs += 1
     
    return ppu, qpu #positive real number

def cap_arr():
    nnode = len(dss.Circuit.AllBusNames())
    Sbase = 1000000.0
    caparr = np.zeros((3, nnode))
    for cap in range(len(dss.Capacitors.AllNames())):
        dss.Capacitors.Name(dss.Capacitors.AllNames()[cap])
        cap_data = dss.CktElement.BusNames()[0].split('.')
        idxbs = dss.Circuit.AllBusNames().index(cap_data[0])
        cap_phases = [0, 0, 0]
        
        for i in range(len(cap_data[1:])):
            cap_phases[int(cap_data[1:][i]) - 1] = 1
        if len(cap_phases) == 0:
            cap_phases = [1, 1, 1]  

        rs = 0
        for ph in range(len(cap_phases)):
            if cap_phases[ph] == 1:
                caparr[ph, idxbs] -= dss.CktElement.Powers()[1::2][rs] * 1e3 / Sbase
                rs += 1

    return caparr #negative real number

def nominal_load_values(t):
    if t == -1:
        var = 1
    else:
        var = (1 + 0.1*np.sin(2*np.pi*0.01*t))
        
    nnode = len(dss.Circuit.AllBusNames())
    Sbase = 1000000.0

    dsskw = np.zeros((3,nnode))
    dsskvar = np.zeros((3,nnode))
    
    for l in dss.Loads.AllNames():
        dss.Loads.Name(l)
        load_phases = [0, 0, 0] #instantiate load phases as all non-existent
        load_data = dss.CktElement.BusNames()[0].split('.')[1:]
        knode = dss.Circuit.AllBusNames().index((dss.CktElement.BusNames()[0].split('.')[0])) #match busname to index
        no_phases = len(load_data)
        for p in load_data:
            phase = int(p)
            load_phases[phase-1] = 1
        if len(load_data) == 0:
            load_phases = [1, 1, 1]        
        realstuff = dss.Loads.kW()
        imagstuff = dss.Loads.kvar()   

        for ph in range(len(load_phases)):      
            if load_phases[ph] == 1:
                dsskw[ph, knode] += realstuff * 1e3 * var / Sbase / no_phases 
                dsskvar[ph, knode] += imagstuff * 1e3 * var / Sbase / no_phases
    
    return dsskw, dsskvar #positive real numbers
 

def nominal_cap_arr():
    nnode = len(dss.Circuit.AllBusNames())
    Sbase = 1000000.0
    caparr = np.zeros((3, nnode))

    for n in range(len(dss.Capacitors.AllNames())):
        dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
        cap_data = dss.CktElement.BusNames()[0].split('.')
        idxbs = dss.Circuit.AllBusNames().index(cap_data[0])
        cap_phases = [0, 0, 0]
        no_phases = len(cap_data[1:])
        if no_phases == 0:
            no_phases = 3
            
        for p in range(len(cap_data[1:])):
            cap_phases[int(cap_data[1:][p]) - 1] = 1
        if len(cap_phases) == 0:
            cap_phases = [1, 1, 1]  
    
        for ph in range(len(cap_phases)):
            if cap_phases[ph] == 1:
                caparr[ph, idxbs] += dss.Capacitors.kvar() * 1e3 / Sbase / no_phases
    return caparr #negative real number

# -------------------------- KCL Functions

def bus_phases():
    dictionary = {}
    for k2 in range(len(dss.Circuit.AllNodeNames())):
        a, b = dss.Circuit.AllNodeNames()[k2].split('.')
        if a in dictionary:
            temp = dictionary[a]
            temp[int(b) - 1] = 1
            dictionary[a] = temp
        elif a not in dictionary:
            dictionary[a] = [0, 0, 0]
            temp = dictionary[a]
            temp[int(b) - 1] = 1
            dictionary[a] = temp
    return dictionary

#line index
def get_line_idx(line):
    return dss.Lines.AllNames().index(line)

#in and out line lists
def linelist(busname):
    in_lines = np.array([])
    out_lines = np.array([])
    for k in range(len(dss.Lines.AllNames())):
        dss.Lines.Name(dss.Lines.AllNames()[k])
        if busname in dss.Lines.Bus1():
            out_lines = np.append(out_lines, dss.Lines.AllNames()[k])
        elif busname in dss.Lines.Bus2():
            in_lines = np.append(in_lines, dss.Lines.AllNames()[k])
    return in_lines,out_lines

# -------------------------- Voltage Regulator and Transformer Parameters

def transformer_regulator_parameters():
    tf_no = len(dss.Transformers.AllNames()) - len(dss.RegControls.AllNames()) #number of transformers
    vr_no = len(dss.RegControls.AllNames()) #number of voltage regulators
  
    tf_bus = np.zeros((5, tf_no), dtype = int) #r1 = in bus, r2 = out bus, r3-5 phases; c = tf
    vr_bus = np.zeros((5, vr_no), dtype = int) 
    gain = np.zeros(vr_no) # gain for voltage regulators
   
    vr_lines = 0
    tf_count = 0
    tf_lines = 0

    for vr in range(vr_no):
        dss.RegControls.Name(dss.RegControls.AllNames()[vr])  
        for i in range(2): 
            #start / end bus
            dss.Transformers.Name(dss.RegControls.Transformer())     
            bus = dss.CktElement.BusNames()[i].split('.')
            vr_bus[i, vr] = int(dss.Circuit.AllBusNames().index(bus[0])) 
        for n in range(len(bus[1:])):                 
            # tap numbers are evenly spaced from -16 to 16 (increments of 0.00625) from 0.9 to 1.1
            # negative sign for voltage ratio purposes                           
            gain_val = (((int(dss.RegControls.TapNumber()) + 16) * (1.1-0.9) / 32) + 0.9)
            gain[vr] = -gain_val
            vr_lines += 1                    
            vr_bus[int(bus[1:][n]) + 1, vr] = int(bus[1:][n]) # phases that exist                   
        if len(bus) == 1: 
            # unspecified phases, assume all 3 exist
            for n in range(1,4): 
                vr_lines += 1
                vr_bus[n+1, vr] = n   
        
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
    
    return tf_bus, vr_bus, tf_lines, vr_lines, tf_count, vr_no, gain

# -------------------------- Slack Bus and KVL Functions

def identify_bus_phases(bus):
        k = np.zeros(3)
        for i in range(1, 4):
            pattern = r"\.%s" % (str(i))
            m = re.findall(pattern, bus)
            if m:
                k[i - 1] = 1
        return k

    
def identify_line_phases(line):
    k = np.zeros(3)
    dss.Lines.Name(line)
    bus = dss.Lines.Bus1()
    for i in range(1, 4):
        pattern = r"\.%s" % (str(i))
        m = re.findall(pattern, bus)
        if m:
            k[i - 1] = 1
    return k


def simple_reg_control(XNR):
    # ptratio = np.array([])
    # tap_winding = np.array([])
    # tap_number = np.array([])
    # forward_vreg = np.array([])
    # kV_reg = np.array([])
    # band = np.array([])
    # reg_bus = np.array([])
    for r in range(len(dss.RegControls.AllNames())):
        dss.RegControls.Name(dss.RegControls.AllNames()[r])
        dss.Transformers.Name(dss.RegControls.AllNames()[r])
        dss.Circuit.SetActiveBus(dss.CktElement.BusNames()[0].split(".")[0])
    
        Vbase = dss.Bus.kVBase() * 1000
        print(Vbase)
        # print(dss.RegControls.AllNames()[r])
        # print(dss.RegControls.PTRatio())
        # print(dss.CktElement.BusNames())
        # print(dss.RegControls.TapWinding())
        # print(dss.RegControls.TapNumber())
        # print(dss.RegControls.ForwardVreg())
        # print(dss.Transformers.kV()) #kv
        # print(dss.RegControls.ForwardBand(), "\n") #band


        # ptratio[r] = dss.RegControls.PTRatio()
        # tap_winding[r] =  dss.RegControls.TapWinding()
        # tap_number[r] = dss.RegControls.TapNumber()
        # forward_vreg[r] = dss.RegControls.ForwardVreg()
        # kV_reg[r] = dss.Transformers.kV()
        # band[r] = dss.RegControls.ForwardBand()
        # reg_bus[r] = dss.CktElement.BusNames()[dss.RegControls.TapWinding() - 1]
        # #need to know the phase and the bus index
        #dss.Transformers.kV() * 1000
        tw = dss.RegControls.TapWinding()
        band = dss.RegControls.ForwardBand()
        fwd_vreg = dss.RegControls.ForwardVreg()

        idxbs = dss.Circuit.AllBusNames().index((dss.CktElement.BusNames()[tw - 1].split('.')[0])) #match busname to index
        ph = dss.CktElement.BusNames()[tw - 1].split('.')[1:] 
        ph_arr = [0, 0, 0]
        if len(ph) == 0:
            ph_arr = [1, 1, 1]
        else:
            for i in ph:
                num = int(i)
                ph_arr[num - 1] = 1

        #NR_voltage =  XNR[2*3*ph + 2*idxbs] / dss.RegControls.PTRatio()
        target_voltage = fwd_vreg
        XNR_final = 0
        for ph in range(len(ph_arr)):
            NR_voltage = np.abs(XNR[2*3*ph + 2*idxbs]) * Vbase
            
            print('NR v', NR_voltage)
            diff = np.abs(NR_voltage - target_voltage)
            print('diff', diff)

            #compare voltage to target voltage
             #assess the difference
            if diff <= band: #converges
                    XNR_final = XNR
            elif diff > band:
                if NR_voltage - target_voltage > 0:
                    if dss.RegControls.TapNumber() >= 16 :
                        print('Tap Number Out of Bounds' )
                        XNR_final = XNR
                 
                    else:
                        print('Increase Tap Number')
                        dss.RegControls.TapNumber(dss.RegControls.TapNumber() + 1)
                        #NR3_timevarying(fn, XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, tol, maxiter, der, capacitance, time_delta, H_reg, G_reg)

                else:
                    if dss.RegControls.TapNumber() <= -16:
                        print('Tap Number Out of Bounds' )
                        XNR_final = XNR
                 
                    else:
                        print('decrease tap number')
                        dss.RegControls.TapNumber(dss.RegControls.TapNumber() - 1)
                        #NR3_timevarying(fn, XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, tol, maxiter, der, capacitance, time_delta, H_reg, G_reg)
        return XNR_final
       
        #make decision about whether to move the tap number or if convergence
        #move the tap number
        #resolve and repeat 


