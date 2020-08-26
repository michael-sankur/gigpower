import numpy as np
import opendssdirect as dss
import re
import sys

class compute_vecmat:

    def __init__(self, XNR, network1, fn, Vslack):
        dss.run_command('Redirect ' + fn)
        dss.Solution.Solve()

        self.all_line_names = dss.Lines.AllNames()
        self.nline = len(self.all_line_names)

        self.all_bus_names = dss.Circuit.AllBusNames()
        self.nnode = len(self.all_bus_names)

        self.all_node_names = dss.Circuit.AllNodeNames()
        self.nn = len(self.all_node_names)

        Vbase = dss.Bus.kVBase() * 1000
        self.Sbase = 1000000.0

        Ibase = self.Sbase/Vbase
        self.Zbase = Vbase/Ibase

        # construct 1 time in line and out line
        self.in_lines = []
        self.out_lines = []

        for k2 in range(1, self.nnode):
            in_lines, out_lines = self.linelist(self.all_bus_names[k2])
            self.in_lines.append(in_lines)
            self.out_lines.append(out_lines)


        #we only need to init this once
        self.bus_phases = self.bus_phases()
        # get all init kW and kvar from OpenDSS, do this only once
        # TODO: create API to set these value manually
        self.buskwkvar = np.zeros((self.nnode, 2, 3))
        for ph in range(3):
            for cplx in range(2):
                for busname in self.all_bus_names:
                    self.buskwkvar[busname, cplx, ph] = self.d_factor(busname, cplx, ph)

        self.line_phases = np.zeros((self.nline,))
        for k, v in enumerate(self.all_line_names):
            self.line_phases[k] = self.identify_line_phases(v)




    def bus_phases(self): #goes through all the buses and saves their phases to a list stored in a dictionary
        #1 if phase exists, 0 o.w.
        #list goes [a, b, c]
        #key is the bus name (without the phase part)
        dictionary = {}
        for k2 in range(self.nn):
            for i in range(1, 4):
                pattern = r"\.%s" % (str(i))

                m = re.findall(pattern, self.all_node_names[k2])
                a, b = self.all_node_names[k2].split('.')
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


    def identify_bus_phases(self, bus): #figures out which phases correspond to the bus
        #returns a list of the r/x matrix places that have those phase/s
        k = np.zeros(3)
        for i in range(1, 4):
            pattern = r"\.%s" % (str(i))
            m = re.findall(pattern, bus)
            if m:
                k[i - 1] = 1
        return k

    def linelist(self, busname): #returns two lists of in and out lines at a bus
        dss.Circuit.SetActiveBus(busname)
        in_lines = np.array([])
        out_lines = np.array([])
        for k in range(self.nline):
            dss.Lines.Name(self.all_line_names[k])
            if busname in dss.Lines.Bus1():
                out_lines = np.append(out_lines, self.all_line_names[k])
            elif busname in dss.Lines.Bus2():
                in_lines = np.append(in_lines, self.all_line_names[k])
        return in_lines, out_lines

    def identify_line_phases(self, line): #figures out which phases correspond to a line
    #(for assigning rmatrix based on line code)
    #returns list of 0's 1's whether or not phase exists in line
        k = np.zeros(3)
        dss.Lines.Name(line)
        bus = dss.Lines.Bus1()
        for i in range(1, 4):
            pattern = r"\.%s" % (str(i))
            m = re.findall(pattern, bus)
            if m:
                k[i - 1] = 1
        return k

    def get_line_idx(self, line): #returns the index of a line as stored in dss.Lines.AllNames()
        k = -1
        for n in range(self.all_line_names):
            if self.all_line_names[n] == line:
                k = n
        return k

    def get_bus_idx(self, bus):
        k = -1
        for n in range(self.all_bus_names): #iterates over all the buses to see which index corresponds to bus
            if self.all_bus_names[n] in bus:
                k = n
        return k

    # this will need to be used only once to get information from openDSS
    def d_factor(self, busname, cplx, ph):
        #factor = np.array([])
        for n in range(len(dss.Loads.AllNames())):
            #dss.Loads.Name(dss.Loads.AllNames()[n])
            if busname in dss.Loads.AllNames()[n]:
                dss.Loads.Name(dss.Loads.AllNames()[n])
                if self.bus_phases[busname][ph] == 0:
                    return 0
                if cplx == 0:
                    return dss.Loads.kW()*1*1e3/self.Sbase
                elif cplx == 1:
                    return  dss.Loads.kvar()*1*1e3/self.Sbase
        return 0