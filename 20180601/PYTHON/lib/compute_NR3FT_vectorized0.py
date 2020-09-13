import numpy as np
def compute_NR3FT_vectorized(X, g_SB, b_SB, H, g, b, nnode):
    FTSUBV = (g_SB @ X) - b_SB

    FTKCL = np.zeros((2*3*(nnode-1), 1))

    for i in range(2*3*(nnode-1)):
        r = (X.T @ (H[i, :, :] @ X)) \
        + (g[i, 0,:] @ X) \
        + b[i, 0,0]
        FTKCL[i, :] = r

    FT = np.r_[FTSUBV,  FTKCL]
    return FT
