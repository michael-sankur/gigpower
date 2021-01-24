import numpy as np
import opendssdirect as dss
import re
import sys
def compute_KCL_matrices(fn, t, der, capacitance, tf_bus, vr_bus, tf_lines, vr_lines):

    dss.run_command('Redirect ' + fn)

    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())
    Sbase = 1000000.0
    
    line_in_idx_vr = range(0, 2*vr_lines, 2)
    line_out_idx_vr = range(1, 2*vr_lines, 2)

    line_idx_tf = range(0, tf_lines)

    # [[3 x nnode array of loads]]
    def load_values():
        load_kw = np.zeros((3, nnode))
        load_kvar = np.zeros((3, nnode))
        if t == -1:
            var = 1
        else:
            var = (1 + 0.1*np.sin(2*np.pi*0.01*t))
        for load in range(len(dss.Loads.AllNames())):
            dss.Loads.Name(dss.Loads.AllNames()[load])
            load_data = dss.CktElement.BusNames()[0].split('.')
            idxbs = dss.Circuit.AllBusNames().index(load_data[0])
            realstuff = dss.CktElement.Powers()[::2]
            imagstuff = dss.CktElement.Powers()[1::2]
            for ph in range(1, len(load_data)):
                load_kw[int(load_data[ph]) - 1, idxbs] += realstuff[ph-1] * 1e3 * var / Sbase 
                load_kvar[int(load_data[ph]) - 1, idxbs] += imagstuff[ph-1] * 1e3 * var/ Sbase
        return load_kw, load_kvar
    load_kw, load_kvar = load_values()

    # [[3 x nnode array of capacitance]]
    def cap_arr():
        caparr = np.zeros((3, nnode))# dtype =complex )
        for n in range(len(dss.Capacitors.AllNames())):
            dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
            cap_data = dss.CktElement.BusNames()[0].split('.')
            idxbs = dss.Circuit.AllBusNames().index(cap_data[0])
            for ph in range(1, len(cap_data)):
                #caparr[int(cap_data[ph]) - 1, idxbs] += dss.Capacitors.kvar() * 1e3 / Sbase / (len(cap_data) - 1)
                caparr[int(cap_data[ph]) - 1, idxbs] += dss.CktElement.Powers()[1::2][0] * 1e3 / Sbase 
        return caparr

    caparr = cap_arr()

    # {bus : [1 x 3 phase existence]}
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

    # ----------Residuals for KCL at a bus (m) ----------
    bp = bus_phases()

    # Zip Parameters
    # Load
    beta_S = 1
    beta_I = 0.0
    beta_Z = 0.0

    # Capacitors
    gamma_S = 0.0
    gamma_I = 0.0
    gamma_Z = 1

    H = np.zeros((2*3*(nnode-1), 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
    g = np.zeros((2*3*(nnode-1), 1, 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
    b = np.zeros((2*3*(nnode-1), 1, 1))

    # Quadratic Terms
    for ph in range(0,3):
        if ph == 0: #set nominal voltage based on phase
            A0 = 1
            B0 = 0
        elif ph == 1:
            A0 = -1/2
            B0 = -1 * np.sqrt(3)/2
        elif ph == 2:
            A0 = -1/2
            B0 = np.sqrt(3)/2
        for k2 in range(1, len(dss.Circuit.AllBusNames())): #skip slack bus
            dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k2]) #set the bus
            in_lines, out_lines = linelist(dss.Circuit.AllBusNames()[k2]) #get in/out lines of bus
            for cplx in range(0,2):
                idxbs = dss.Circuit.AllBusNames().index(dss.Circuit.AllBusNames()[0])
                if cplx == 0:
                    load_val = load_kw[ph][idxbs]
                    cap_val = 0
                else:
                    load_val = load_kvar[ph][idxbs]
                    cap_val = caparr[ph][idxbs]
                gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))]) #some derivatives
                hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
                                    [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])
                available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                if available_phases[ph] == 1:                 #quadratic terms
                    H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2] = \
                                                                                                    -load_val * (beta_Z + (0.5 * beta_I* hessian_mag[0][0])) + \
                                                                                                    cap_val * (gamma_Z + (gamma_I * hessian_mag[0][0]))# TE replace assignment w/ -load_val * beta_Z; #a**2
                    H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2 + 1] = \
                                                                                                    -load_val * (beta_Z + (0.5 * beta_I * hessian_mag[1][1])) + \
                                                                                                    cap_val * (gamma_Z + (gamma_I * hessian_mag[1][1]))# TE replace assignment w/ -load_val * beta_Z; #b**2
                        #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1] = -load_val * beta_I * hessian_mag[0][1] / 2 #remove for TE
                        #H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2] =  -load_val * beta_I * hessian_mag[1][0] / 2 #remove for TE

                for i in range(len(in_lines)): #fill in H for the inlines            
                    line_idx = get_line_idx(in_lines[i])
                    if available_phases[ph] == 1:
                        if cplx == 0: #real residual
                            #A_m and C_lm
                            H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx] = 1/2
                            H[2*ph*(nnode-1) + (k2-1)*2 + cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2] = 1/2
                            #B_m and D_lm
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = 1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = 1/2
                        if cplx == 1: #imaginary residual
                            # #A_m, D_lm
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = -1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = -1/2
                            #B_m and C_lm
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx] = 1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = 1/2
            
                for j in range(len(out_lines)): #fill in H for the outlines    
                    line_idx = get_line_idx(out_lines[j])
                    if available_phases[ph] == 1:
                        if cplx == 0: #real residual
                            #A_m and C_mn
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx] = -1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2] = -1/2
                            #B_m and D_mn
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1] = -1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = -1/2
                        if cplx == 1: #imaginary residual
                            #A_m and D_mn
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1]= 1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+ cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = 1/2
                            #C_m and B_mn
                            H[2*ph*(nnode-1) + (k2-1)*2+cplx][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx] = -1/2
                            H[2*ph*(nnode-1) + (k2-1)*2+cplx][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = -1/2

    # Transformer KCL
    tf_no = len(dss.Transformers.AllNames()) - len(dss.RegControls.AllNames())
    count_tf = 0
    count_tf2 = 0
    for i in range(tf_no):
        for ph in range(0,3):     
            k2 = int(tf_bus[1, i]) #in bus index of transformer [out bus: a0] --line-a0-a1-- [in bus: a1]
            if k2 == 0: #if source bus, need to skip line
                count_tf += 1
            if k2 != 0 and tf_bus[ph + 2, i] != 0: #if not source bus, perform KCL             
                line_idx = line_idx_tf[count_tf]                  
                #A_m and C_lm
                H[2*ph*(nnode-1) + (k2-1)*2 ][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*line_idx] = 1/2
                H[2*ph*(nnode-1) + (k2-1)*2 ][2*3*(nnode+nline) + 2*line_idx][2*(nnode)*ph + 2*k2] = 1/2
                #B_m and D_lm
                H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*line_idx + 1] = 1/2
                H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = 1/2
                        
                #A_m, D_lm
                H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*line_idx + 1] = -1/2
                H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*3*(nnode+nline) + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = -1/2
                #B_m and C_lm
                H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*line_idx] = 1/2
                H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*3*(nnode+nline) + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = 1/2
                count_tf += 1 #go to next line
        
    for j in range(tf_no): #fill in H for the outlines  
        for ph in range(0,3):                 
            k2 = int(tf_bus[0, j]) #out bus index of transformer
            if k2 == 0: 
                count_tf2 += 1
            if k2 != 0 and tf_bus[ph + 2, j] != 0:
                line_idx = line_idx_tf[count_tf2]                    
                #real residual
                #A_m and C_mn        
                H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*line_idx] = -1/2
                H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*line_idx][2*(nnode)*ph + 2*k2] = -1/2
                #B_m and D_mn
                H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*line_idx + 1] = -1/2
                H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = -1/2                       
                
                #imaginary residual               
                #A_m and D_mn
                H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*line_idx + 1]= 1/2
                H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*3*(nnode+nline) + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = 1/2
                #C_m and B_mn
                H[2*ph*(nnode-1) + (k2-1)*2+1][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*line_idx] = -1/2
                H[2*ph*(nnode-1) + (k2-1)*2+1][2*3*(nnode+nline) + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = -1/2
                count_tf2+=1

    #Voltage Regulator KCL 
    vr_count = len(dss.RegControls.AllNames())
    count_vr = 0
    count_vr2 = 0   
    if vr_count > 0:
        for i in range(vr_count): #in lines 
            for ph in range(0,3):                    
                k2 = int(vr_bus[1, i])
                if k2 == 0: 
                    count_vr += 1
                if k2 != 0 and vr_bus[ph + 2, i] != 0:    
                    line_idx = line_in_idx_vr[count_vr]   
                    #real residual
                    #A_m and C_lm
                    H[2*ph*(nnode-1) + (k2-1)*2 + 0][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx] = 1/2
                    H[2*ph*(nnode-1) + (k2-1)*2 + 0][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx][2*(nnode)*ph + 2*k2] = 1/2
                    #B_m and D_lm
                    H[2*ph*(nnode-1) + (k2-1)*2+ 0][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1] = 1/2
                    H[2*ph*(nnode-1) + (k2-1)*2+ 0][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = 1/2
                    #imaginary residual                        
                    # #A_m, D_lm
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1] = -1/2
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = -1/2
                    #B_m and C_lm
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx] = 1/2
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = 1/2
                    count_vr += 1

        for j in range(vr_count): #fill in H for the outlines       
            for ph in range(0,3):                 
                k2 = int(vr_bus[0, j])
                if k2 == 0: 
                    count_vr2 += 1
                if k2 != 0 and vr_bus[ph + 2, j] != 0:
                    line_idx = line_out_idx_vr[count_vr2]     
                    #real residual
                    #A_m and C_mn
                    H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx] = -1/2
                    H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx][2*(nnode)*ph + 2*k2] = -1/2
                    #B_m and D_mn
                    H[2*ph*(nnode-1) + (k2-1)*2][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1] = -1/2
                    H[2*ph*(nnode-1) + (k2-1)*2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1] = -1/2
                    #imaginary residual
                    #A_m and D_mn
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*(nnode)*ph + 2*k2][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1] = 1/2
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx + 1][2*(nnode)*ph + 2*k2] = 1/2
                    #C_m and B_mn
                    H[2*ph*(nnode-1) + (k2-1)*2+ 1][2*(nnode)*ph + 2*k2 + 1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx] = -1/2
                    H[2*ph*(nnode-1) + (k2-1)*2 +1][2*3*(nnode+nline) + 2*tf_lines + 2*line_idx][2*(nnode)*ph + 2*k2 + 1] = -1/2
                    count_vr2 += 1

    
    #Linear Term & Constant Term
    for ph in range(0,3):
        if ph == 0: #set nominal voltage based on phase
            A0 = 1
            B0 = 0
        elif ph == 1:
            A0 = -1/2
            B0 = -1 * np.sqrt(3)/2
        elif ph == 2:
            A0 = -1/2
            B0 = np.sqrt(3)/2
        for k2 in range(1, len(dss.Circuit.AllBusNames())):
            for cplx in range(0,2):
                available_phases = bp[dss.Circuit.AllBusNames()[k2]] #phase array at specific bus
                idxbs = dss.Circuit.AllBusNames().index(dss.Circuit.AllBusNames()[k2])
                if cplx == 0:
                    load_val = load_kw[ph][idxbs]
                else:
                    load_val = load_kvar[ph][idxbs]

                # linear terms
                g_temp = np.zeros(2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines) #preallocate g
                if available_phases[ph] == 0: #if phase does not exist
                    g_temp[2*(ph)*nnode + 2*k2 + cplx] = 1
                g[2*(nnode-1)*ph + 2*(k2-1) + cplx, 0,:] = g_temp

                # Constant terms
                if cplx == 0:
                    if der.real != 0:
                        b_factor = der.real
                        b_factor = 0
                    else:
                        b_factor = 0
                elif cplx == 1:
                    if capacitance != 0 or der.imag != 0:
                        b_factor = capacitance - der.imag
                        b_factor = 0
                    else:                    
                        cap_val = caparr[ph][k2]                      
                else:
                    b_factor = 0

                if available_phases[ph] == 0: #if phase does not exist at bus, set b = 0
                    b_temp = 0
                else:
                    b_temp = (-load_val * beta_S) + (cap_val * gamma_S) #TE version
                   
                    # b_temp = -load_val * (beta_S \
                    # + (beta_I) * (((hessian_mag[0][1] * A0 * B0) + ((1/2)*hessian_mag[0][0] * ((A0)**2)) + ((1/2)*hessian_mag[1][1] * (B0**2))) \
                    # -  (A0 * gradient_mag[0] + B0* gradient_mag[1]) \
                    # +  (A0**2 + B0**2) ** (1/2))) \
                    # + b_factor #calculate out the constant term in the residual

                b[2*(nnode-1)*ph + 2*(k2-1) + cplx][0][0] = b_temp #store the in the b matrix


    return H, g, b
