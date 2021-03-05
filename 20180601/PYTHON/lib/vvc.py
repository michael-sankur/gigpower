import opendssdirect as dss
class voltVARControl:
    def __init__(self, voltage_points, minQ, maxQ, busName, phase):
        # voltage breakpoints
        self.V1 = voltage_points[0]
        self.V2 = voltage_points[1]
        self.V3 = voltage_points[2]
        self.V4 = voltage_points[3]
        self.voltage_points = voltage_points
        self.minQ = minQ
        self.maxQ = maxQ
        self.prevQ_list = []
        self.prevQ = 0
        self.busName = busName
        self.phase = phase

    def __str__(self):
        return self.busName

    def get_Q(self, V_pu):
        if V_pu <= self.V1:
            self.prevQ = self.maxQ
        elif V_pu <= self.V2:
            slope =  (- self.maxQ) / (self.V2 - self.V1)
            b = - (slope * self.V2)
            self.prevQ = slope * V_pu + b
        elif V_pu <= self.V3:
            self.prevQ = 0
        elif V_pu <= self.V4:
            slope = (self.minQ) / (self.V4- self.V3)
            b = - (slope * self.V3)
            self.prevQ = slope * V_pu + b
        else:
            self.prevQ = self.minQ
        
        self.prevQ_list.append(self.prevQ)
        return self.prevQ
    
    def get_prevQ(self):
        return self.prevQ

    def get_busName(self):
        return self.busName
    
    def get_phase(self):
        return self.phase