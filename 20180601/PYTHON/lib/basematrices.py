import numpy as np
from lib.compute_vecmat import compute_vecmat
from lib.compute_KCL_matrices import compute_KCL_matrices
import opendssdirect as dss
import re
import time
def basematrices(fn, slacknode, Vslack, V0, I0):

    dss.run_command('Redirect ' + fn)
      
    tf_no = len(dss.Transformers.AllNames()) #number of transformers
    vr_no = len(dss.RegControls.AllNames()) #number of voltage regulators
  
    tf_bus = np.zeros((5, tf_no), dtype = int) #r1 = in bus, r2 = out bus, r3-5 phases; c = tf
    vr_bus = np.zeros((5, vr_no), dtype = int) 
 
    vr_count = 0 
    vr_lines = 0
    tf_count = 0
    tf_lines = 0

    for vr in range(len(dss.RegControls.AllNames())):
        dss.RegControls.Name(dss.RegControls.AllNames()[vr])  
        for i in range(2):
            dss.Transformers.Name(dss.RegControls.Transformer())     
            bus = dss.CktElement.BusNames()[i].split('.')
            vr_bus[i, vr_count] = int(dss.Circuit.AllBusNames().index(bus[0])) #start / end bus
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
            continue
        tf_bus[0:2, tf_count] = tf_bus_temp[:, 0]
        for n in range(len(bus[1:])):                 
            tf_lines += 1                             
            tf_bus[int(bus[1:][n]) + 1, tf_count] = int(bus[1:][n])   
        if len(bus) == 1:
            for k in range(1,4):
                tf_lines += 1
                tf_bus[k+1, tf_count] = k             
        tf_count += 1

    idx = np.argwhere(np.all(tf_bus[:, :] == 0, axis=0))
    tf_bus = np.delete(tf_bus, idx, axis=1)
    
    nline = len(dss.Lines.AllNames())  
    nnode = len(dss.Circuit.AllBusNames())
    XNR = np.zeros((2*3*(nnode + nline) + 2*tf_lines + 2*2*vr_lines,1))


    # intialize node voltage portion of XNR
    if V0 == None or len(V0) == 0:
        for ph in range(0,3):
            for k1 in range(0,nnode):
                XNR[2*ph*nnode + 2*k1] = Vslack[ph].real
                XNR[2*ph*nnode + 2*k1+1] = Vslack[ph].imag

    # If initial V is given (usually from CVX)
    elif len(V0) != 0:
        for ph in range(0,3):
            for k1 in range(0,nnode):
                XNR[2*ph*nnode + 2*k1] = V0[ph,k1].real
                XNR[2*ph*nnode + 2*k1+1] = V0[ph,k1].imag

    # intialize line current portion of XNR
    if I0 == None or len(I0) == 0:
        XNR[(2*3*nnode):] = 0.0*np.ones((6*nline + 2*tf_lines + 2*2*vr_lines ,1))

    # If initial I is given
    elif len(I0) != 0:
        for ph in range(0,3):
            for k1 in range(0,len(dss.Lines.AllNames())):
                XNR[(2*3*nnode) + 2*ph*nline + 2*k1] = I0[ph,k1].real
                XNR[(2*3*nnode) + 2*ph*nline + 2*k1+1] = I0[ph,k1].imag
        XNR[(2*3*nnode + 2*3*nline):] = np.zeros((len(XNR) - 2*3*nnode - 2*3*nline), 1) 
    
    # generate static matrices
    XNR, g_SB, b_SB, G_KVL, b_KVL, H_reg, G_reg = compute_vecmat(XNR, fn, Vslack)
    # generate non-static matrices
    H, g, b = compute_KCL_matrices(fn, -1, 0, 0)
    return XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg
