import opendssdirect as dss
from lindist3flow import LinDist3Flow
import numpy as np
import time
dss_path = '/home/toanngo/Documents/GitHub/LinDist3Flow/lindist3flow/data/compare_opendss_05node_threephase_unbalanced_oscillation_03.dss'
#dss_path = '/home/toanngo/Documents/GitHub/LinDist3Flow/lindist3flow/data/IEEE_13Node_Modified_01.dss'
#dss_path = '/home/toanngo/Documents/GitHub/LinDist3Flow/lindist3flow/data/IEEE_13Node_Modified_01.dss'
#dss_path = '/home/toanngo/Documents/GitHub/LinDist3Flow/lindist3flow/data/06node_threephase_unbalanced.dss'
dss.run_command('Redirect ' + dss_path)
nnode = len(dss.Circuit.AllBusNames())
nline = len(dss.Lines.AllNames())


# loads
# Store the loads
slackidx = 0
Vslack = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])
times = np.linspace(0, 2*np.pi, 5)

kW_list = np.array([])
kvar_list = np.array([])
for k in range(len(dss.Loads.AllNames())):
    dss.Loads.Name(dss.Loads.AllNames()[k])
    kW_list = np.append(kW_list, dss.Loads.kW())
    kvar_list = np.append(kvar_list, dss.Loads.kvar())
all_load_names = dss.Loads.AllNames()

VNR01 = np.zeros((len(times), 3, nnode), dtype = "complex")
INR01 =  np.zeros((len(times), 3, nline), dtype = "complex")
STXNR01 = np.zeros((len(times), 3, nline), dtype = "complex")
SRXNR01 = np.zeros((len(times), 3, nline), dtype = "complex")
iNR01 = np.zeros((len(times), 3, nnode), dtype = "complex")
sNR01 = np.zeros((len(times), 3, nnode), dtype = "complex")

LinDist3Flow = LinDist3Flow(dss_path, slackidx, Vslack)
start_time = time.time()
for i in range(len(times)-1):

    VNR, INR, iNR, sNR, STXNR, SRXNR = LinDist3Flow.solve()
    # Time-varying load
    for k,v in enumerate(all_load_names):
        LinDist3Flow.set_load_kw(v, kW_list[k]* (1 + 0.1*np.sin(2*np.pi*0.01*times[i+1])))
        LinDist3Flow.set_load_kvar(v, kvar_list[k] * (1 + 0.1*np.sin(2*np.pi*0.01*times[i+1])))
    #VNR01[i, :, :] = np.reshape(VNR, (3, nnode))
    #INR01[i, :, :] = np.reshape(INR, (3, nline))
    #STXNR01[i, :, :] = np.reshape(STXNR, (3, nline))
    #SRXNR01[i, :, :] = np.reshape(SRXNR, (3, nline))
    #iNR01[i, :, :] = np.reshape(iNR, (3, nnode))
    #sNR01[i, :, :] = np.reshape(sNR, (3, nnode))
    print(VNR)
print(time.time() - start_time)