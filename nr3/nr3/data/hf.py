import opendssdirect as dss
import time
import re
import numpy as np
def bus_phases(): #goes through all the buses and saves their phases to a list stored in a dictionary
#1 if phase exists, 0 o.w.
#list goes [a, b, c]
#key is the bus name (without the phase part)
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
    return dictionary

def cap_arr():
    caparr = np.zeros((3, nnode))
    for n in range(len(dss.Capacitors.AllNames())):
        dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
        cap_data = dss.CktElement.BusNames()[0].split('.')

        idxbs = dss.Circuit.AllBusNames().index(cap_data[0])
        for ph in range(1, len(cap_data)):
            caparr[int(cap_data[ph]) - 1, idxbs] += dss.Capacitors.kvar() * 1e3 / Sbase / (len(cap_data) - 1)
    return caparr

def load_order_f():
    load_order = {}
    for n in range(len(dss.Loads.AllNames())):
        dss.Loads.Name(dss.Loads.AllNames()[n])
        pattern =  r"(\w+)\."
        load_bus = re.findall(pattern, dss.CktElement.BusNames()[0])
        if load_bus[0] not in load_order:
            load_order[load_bus[0]] = 1
        elif load_bus[0] in load_order:
            load_order[load_bus[0]] += 1
    return load_order
load_order_list = load_order_f()

def load_values():
    load_ph_arr = np.zeros((nnode, max(load_order_list.values()), 3))
    load_kw_arr_ph = np.zeros((3, nnode))
    load_kvar_arr_ph = np.zeros((3, nnode))
    if t == -1:
        var = 1
    else:
        var = (1 + 0.1*np.sin(2*np.pi*0.01*t))
    for load in range(len(dss.Loads.AllNames())):
        dss.Loads.Name(dss.Loads.AllNames()[load])
        pattern =  r"(\w+)\."
        load_bus = re.findall(pattern, dss.CktElement.BusNames()[0])
        load_ph_arr_temp = [0, 0, 0]
        for i in range(1, 4):
            pattern = r"\.%s" % (str(i))
            load_ph = re.findall(pattern, dss.CktElement.BusNames()[0])
            if load_ph:
                load_ph_arr_temp[i - 1] = 1
        for j in range(max(load_order_list.values())):
            idxbs = dss.Circuit.AllBusNames().index(load_bus[0])
            if np.all(load_ph_arr[idxbs, j,:] == [0, 0, 0]):
                load_ph_arr[idxbs, j, :] = load_ph_arr_temp
                for i in range(len(load_ph_arr_temp)):
                    if load_ph_arr_temp[i] == 1:
                        load_kw_arr_ph[i,idxbs] += dss.Loads.kW()*1e3*var / Sbase  / sum(load_ph_arr_temp)
                        load_kvar_arr_ph[i,idxbs] += dss.Loads.kvar()*1e3*var / Sbase / sum(load_ph_arr_temp)
                break
    return load_kw_arr_ph, load_kvar_arr_ph
