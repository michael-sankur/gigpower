import numpy as np
def compute_NR3JT_vectorized(X, g_SB, G_KVL, H, g, nnode, nline):
    JSUBV = g_SB
    JKVL = G_KVL
    JKCL = np.zeros((2*3*(nnode-1), 2*3*(nnode+nline)))

    for i in range(2*3*(nnode-1)):
        r = (2 * (X.T @ H[:, :, i])) \
        + (g[0,:,i])
        JKCL[i,:] = r

    JT = np.r_[JSUBV, JKVL, JKCL]

    return JT
