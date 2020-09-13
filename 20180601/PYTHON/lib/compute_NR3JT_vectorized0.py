import numpy as np
def compute_NR3JT_vectorized(X, g_SB, H, g, nnode, nline):
    JSUBV = g_SB
    JKCL = np.zeros((2*3*(nnode-1), 2*3*(nnode)))

    for i in range(2*3*(nnode-1)):
        r = (2 * (X.T @ H[i, :, :])) \
        + (g[i, 0, :])
        JKCL[i,:] = r


    JT = np.r_[JSUBV,  JKCL]

    return JT
