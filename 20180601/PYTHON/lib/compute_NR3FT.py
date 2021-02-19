import numpy as np
def compute_NR3FT(X, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode, nline, H_reg, G_reg, vr_lines):
    FTSUBV = (g_SB @ X) - b_SB
    FTKVL = (G_KVL @ X) - b_KVL

    FTKCL = np.zeros((2*3*(nnode-1), 1))
    for i in range(2*3*(nnode-1)):
        r = (X.T @ (H[i, :, :] @ X)) \
        + (g[i, 0,:] @ X) \
        + b[i, 0,0]
        FTKCL[i, :] = r
   
    FTVR = np.zeros((2*vr_lines, 1)) #need to fix 
    for i in range(2*vr_lines):
        r = X.T @ (H_reg[i, :, :] @ X)       
        FTVR[i, :] = r

    FTVR2 = G_reg @ X
   
    FT = np.r_[FTSUBV, FTKVL, FTKCL, FTVR, FTVR2]

    return FT
