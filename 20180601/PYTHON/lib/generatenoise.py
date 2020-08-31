import numpy as np
import random

def generate_spu_noise():
    # n - spu random noise
    n = np.random.normal(0, 1)
    with open("spu_noise.txt", "a") as txt_file:
        txt_file.write(str(n) + "\n")
    return n

def generate_pQ_z_I_noise(ts, nnode):
    # Generate the noise & save to text files
    
    pQ_noise = np.zeros((ts, 3, nnode))
    Z_noise = np.zeros((ts, 3, nnode))
    I_noise = np.zeros((ts, 3, nnode))
       
    for kt in range(ts):
        for i in range(3):
            for j in range(nnode):
                # generate noise
                pQ = int(random.random() * (1/2) * 1000)/1000
                I = int(random.random() * (1/2)*1000)/1000
                Z = 1 - pQ - I

                # store noise to compare w matlab output
                pQ_noise[kt][i][j] = pQ
                Z_noise[kt][i][j] = Z
                I_noise[kt][i][j] = I
    
    # write noise to text files
    with open("aPQ.txt", "w") as txt_file:
        for num in pQ_noise:
            for s in num:
                for j in s:
                    txt_file.write(str(j) + " ")
                txt_file.write("\n")
    
    with open("aI.txt", "w") as txt_file:
        for num in I_noise:
            for s in num:
                for j in s:
                    txt_file.write(str(j) + " ")
                txt_file.write("\n")
    with open("aZ.txt", "w") as txt_file:
        for num in Z_noise:
            for s in num:
                for j in s:
                    txt_file.write(str(j) + " ")
                txt_file.write("\n")
    return pQ_noise, Z_noise, I_noise

def read_a_noises(nnode, ts):
    with open(r'C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\aI.txt', 'r') as f:
        l = [[num for num in line.split(' ')] for line in f]
    l = np.array(l)
    aI = np.delete(l, -1, axis = 1).reshape(ts, 3, nnode)
    
    with open(r'C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\aZ.txt', 'r') as f:
        l = [[num for num in line.split(' ')] for line in f]
    l = np.array(l)
    aZ = np.delete(l, -1, axis = 1).reshape(ts, 3, nnode)

    with open(r'C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\aPQ.txt', 'r') as f:
        l = [[num for num in line.split(' ')] for line in f]
    l = np.array(l)
    aPQ = np.delete(l, -1, axis = 1).reshape(ts, 3, nnode)
    
    return aPQ, aZ, aI

def add_noise(network1, kt, nnode, pQ_noise,  Z_noise, I_noise):
    # network1 - network being adjusted
    # kt       - timestep (ts) of corresponding network
    # nnode    - number of nodes
    # pQ_noise - ts x phases x nnode matrix of noise
    # Z_noise  - ts x phases x nnode matrix of noise
    # I_noise  - ts x phases x nnode matrix of noise

    for i in range(3):
        for j in range(nnode):
            # adjust aPQ, aI, and aZ
            network1.loads.aPQ[i][j] = pQ_noise[kt][i][j]*(network1.loads.spu[i][j] != 0)
            network1.loads.aI[i][j] = I_noise[kt][i][j]*(network1.loads.spu[i][j] != 0)
            network1.loads.aZ[i][j] = Z_noise[kt][i][j]*(network1.loads.spu[i][j] != 0)

    return network1
    return pQ_noise, Z_noise, I_noise