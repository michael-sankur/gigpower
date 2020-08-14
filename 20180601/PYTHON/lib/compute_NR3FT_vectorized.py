import numpy as np
def compute_NR3FT_vectorized(X, g_SB, b_SB, G_KVL, b_KVL, H, g, b, nnode):
    FTSUBV = (g_SB @ X) - b_SB
    FTKVL = (G_KVL @ X) - b_KVL

    FTKCL = np.zeros((2*3*(nnode-1), 1))

    for i in range(2*3*(nnode-1)):
        r = (X.T @ (H[:, :, i] @ X)) \
        + (g[0,:,i] @ X) \
        + b[0,0,i]
        FTKCL[i,:] = r

    FT = np.r_[FTSUBV, FTKVL, FTKCL]
    a_file = open("vectorized.txt", "a+")
    a_file.write('FTSUBV: \n')
    for row in FTSUBV:
        np.savetxt(a_file, row)
    a_file.write('FTKVL: \n')
    for row in FTKVL:
        np.savetxt(a_file, row)
    a_file.write('FTKCL: \n')
    for row in FTKCL:
        np.savetxt(a_file, row)

    a_file.close()
    return FT
