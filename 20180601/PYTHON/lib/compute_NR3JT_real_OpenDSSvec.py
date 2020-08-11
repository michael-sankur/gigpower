import numpy as np
import opendssdirect as dss

def compute_NR3JT_real_function(g_SB, G_KVL, H, X_T, g):

    dss.run_command('Redirect compare_opendss_05node_threephase_unbalanced_oscillation_03.dss')
    dss.Solution.Solve()
    nline = len(dss.Lines.AllNames())
    nnode = len(dss.Circuit.AllBusNames())

    # Jacobian for slack node voltage
    JSUBV = g_SB
    # Jacobian for KVL across lines (m,n)
    JKVL = G_KVL
    # Jacobian for KCL at node m
    JKCL = np.zeros((2*3*(nnode-1), 2*3*(nnode+nline)))

    for i in range(2*3*(nnode-1)):
        r = (2 * (X_T[:, :, i] @ H[:, :, i])) \
        + (g[0,:,i])
        JKCL[i,:] = r

    JT = np.r_[JSUBV, JKVL, JKCL]

    return JT
