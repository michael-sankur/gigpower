import numpy as np
def compute_NR3JT_vectorized(X, g_SB, G_KVL, H, g, nnode, nline):
    JSUBV = g_SB
    JKVL = G_KVL
    JKCL = np.zeros((2*3*(nnode-1), 2*3*(nnode)))
     #np.zeros((2*3*(nnode-1),2*3*(nnode + nline)))

    for i in range(2*3*(nnode-1)):
        r = (2 * (X.T @ H[i, :, :])) \
        + (g[i, 0, :])
        JKCL[i,:] = r


    JT = np.r_[JSUBV, JKVL, JKCL]

    # a_file = open("vectorized.txt", "a+")
    # a_file.write('JSUBV: \n')
    # for row in JSUBV:
    #     np.savetxt(a_file, row)
    # a_file.write('JKVL: \n ')
    # for row in JKVL:
    #     np.savetxt(a_file, row)
    # a_file.write('JKCL: \n')
    # for row in JKCL:
    #     np.savetxt(a_file, row)
    #
    # a_file.close()
    # print('jsubv')
    # print(JSUBV.shape)
    # print('jkvl')
    # print(JKVL.shape)
    # print('jkcl: ')
    # print(JKCL.shape)

    return JT
