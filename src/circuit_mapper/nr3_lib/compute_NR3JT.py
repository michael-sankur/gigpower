import numpy as np
def compute_NR3JT(X, g_SB, G_KVL, H, g, nnode, nline, H_reg, G_reg, tf_lines, vr_lines):
    JSUBV = g_SB
    JKVL = G_KVL

    JKCL = np.zeros((2*3*(nnode-1), 2*3*(nnode+nline) + 2*tf_lines + 2*2*vr_lines))
    for i in range(2*3*(nnode-1)):
        r = (2 * (X.T @ H[i, :, :])) \
        + (g[i, 0, :])
        JKCL[i,:] = r

    JVR = np.zeros((2*vr_lines, 2*3*(nnode+nline)+ 2*tf_lines + 2*2*vr_lines))
    for i in range(2*vr_lines):
        r = (2 * (X.T @ H_reg[i, :, :]))      
        JVR[i, :] = r

    JVR2 = G_reg

    JT = np.r_[JSUBV, JKVL, JKCL, JVR, JVR2]

    return JT
