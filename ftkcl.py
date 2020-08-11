
H = np.zeros((2 * 3 * (nnode + nline), 2 * 3* (nnode + nline), 2*3*(nnode-1)))
g = np.zeros((1, 2*3*(nnode+nline), 2*3*(nnode-1)))
b = np.zeros((1, 1, 2*3*(nnode-1)))

for ph in range(0,3):
    if ph == 0: #set nominal voltage based on phase
        A0 = 1
        B0 = 0
    elif ph == 1:
        A0 = -1/2
        B0 = -1 * np.sqrt(3)/2
    elif ph == 2:
        A0 = -1/2
        B0 = np.sqrt(3)/2
    for k2 in range(1, len(dss.Circuit.AllBusNames())): #skip slack bus
        if k2 == 0: #skip slackidx
            continue
        dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k2])
        in_lines, out_lines = linelist(dss.Circuit.AllBusNames()[k2]) #get in/out lines of bus
        for cplx in range(0,2):
            load_val = d_factor(dss.Circuit.AllBusNames()[k2], cplx)
            gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))]) #some derivatives
            hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
                                [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])

            #quadratic terms
            H[2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2][2*ph*(nnode-1) + (k2-1)*2] = -load_val * (beta_Z + (0.5 * beta_I* hessian_mag[0][0])) #a**2
            H[2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2 + 1][2*ph*(nnode-1) + (k2-1)*2] = -load_val * (beta_Z + (0.5 * beta_I * hessian_mag[1][1])) #b**2
            H[2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1][2*ph*(nnode-1) + (k2-1)*2] = -load_val * beta_I * hessian_mag[0][1] #cross quad. terms in taylor exp
            H[2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2][2*ph*(nnode-1) + (k2-1)*2] =  -load_val * beta_I * hessian_mag[0][1]

            for i in range(len(in_lines)): #fill in H for the inlines
                dss.Lines.Name(in_lines[i])
                line_idx = get_line_idx(in_lines[i])
                if cplx == 0: #real residual
                    #A_m and C_lm
                    H[2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*ph*(nnode-1) + (k2-1)*2] = 1/2
                    H[2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2][2*ph*(nnode-1) + (k2-1)*2] = 1/2
                    #B_m and D_lm
                    H[2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*ph*(nnode-1) + (k2-1)*2] = 1/2
                    H[2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1][2*ph*(nnode-1) + (k2-1)*2] = 1/2
                if cplx == 1: #complex residual
                    #A_m, D_lm
                    H[2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*ph*(nnode-1) + (k2-1)*2] = -1/2
                    H[2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2][2*ph*(nnode-1) + (k2-1)*2] = -1/2
                    #B_m and C_lm
                    H[2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*ph*(nnode-1) + (k2-1)*2] = 1/2
                    H[2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1][2*ph*(nnode-1) + (k2-1)*2] = 1/2

            for j in range(len(out_lines)): #fill in H for the outlines
                dss.Lines.Name(out_lines[j])
                line_idx = get_line_idx(out_lines[j])
                k = get_bus_idx(dss.Lines.Bus2())
                if cplx == 0:
                    #A_m and C_mn
                    H[2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*ph*(nnode-1) + (k2-1)*2] = -1/2
                    H[2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2][2*ph*(nnode-1) + (k2-1)*2] = -1/2
                    #B_m and D_mn
                    H[2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*ph*(nnode-1) + (k2-1)*2] = -1/2
                    H[2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1][2*ph*(nnode-1) + (k2-1)*2] = -1/2
                if cplx == 1:
                    #A_m and D_mn
                    H[2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*ph*(nnode-1) + (k2-1)*2] = 1/2
                    H[2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2][2*ph*(nnode-1) + (k2-1)*2] = 1/2
                    #C_m and B_mn
                    H[2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*ph*(nnode-1) + (k2-1)*2] = -1/2
                    H[2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1][2*ph*(nnode-1) + (k2-1)*2] = -1/2

#Linear Term
for ph in range(0,3):
    if ph == 0: #set nominal voltage based on phase
        A0 = 1
        B0 = 0
    elif ph == 1:
        A0 = -1/2
        B0 = -1 * np.sqrt(3)/2
    elif ph == 2:
        A0 = -1/2
        B0 = np.sqrt(3)/2
    # for k3 in range(len(dss.Lines.AllNames())):
    for k2 in range(len(dss.Circuit.AllBusNames())):
        if k2 == 0: #skip slackidx
            continue
        for cplx in range(0,2):
            load_val = d_factor(dss.Circuit.AllBusNames(), cplx)
            #linear terms
            g_temp = np.zeros(len(X))
            g_temp[2*ph*nnode+ 2 * k2] = -load_val* beta_I* ((1/2 * (-2 * A0 * hessian_mag[0][0] - 2 * B0 * hessian_mag[0][1])) \
                                   +   gradient_mag[0])
            g_temp[2*ph*nnode+ 2 * k2 + 1] = -load_val * beta_I * ((1/2 * (-2* A0 *hessian_mag[0][1] - 2 * B0 * hessian_mag[1][1])) \
                                       +  gradient_mag[1])
            g[0,:,2*(nnode-1)*ph + 2*(k2-1) + cplx] = g_temp

#Constant Term
for ph in range(0,3):
    if ph == 0: #set nominal voltage based on phase
        A0 = 1
        B0 = 0
    elif ph == 1:
        A0 = -1/2
        B0 = -1 * np.sqrt(3)/2
    elif ph == 2:
        A0 = -1/2
        B0 = np.sqrt(3)/2

    for k2 in range(len(dss.Circuit.AllBusNames())):
        if k2 == 0:
            continue
        for cplx in range(0,2):
            load_val = d_factor(dss.Circuit.AllBusNames()[k2], cplx)
            #constant terms
            b_factor = 0
            Sk = dss.CktElement.Powers() #retrieve powers
            if cplx == 1:
                b_factor = dss.Capacitors.kvar() - Sk[1] #depends on if it's real or im residual
            elif cplx == 0:
                b_factor = - Sk[0]
            b_temp = -load_val * (beta_S \
                + (beta_I) * (hessian_mag[0][1] * A0 * B0 + (1/2)*hessian_mag[0][0] * A0**2 + (1/2)*hessian_mag[1][1] * B0**2) \
                - beta_I * (A0 * gradient_mag[0] +B0* gradient_mag[1]) \
                + beta_I * (A0**2 + B0**2) ** (1/2)) \
                + b_factor #calculate out the constant term in the residual
            b[0][0][2*(nnode-1)*ph + 2*(k2-1) + cplx] = b_temp #store the in the b matrix


Y = X.reshape(-1, 1) #(72,1)
enlarged_X = np.zeros((2*3*(nline+nnode), 1, 2*3*(nnode-1)))
X_T= np.zeros((1, 2*3*(nline+nnode), 2*3*(nnode-1)))
for n in range(2*3*nline):
    enlarged_X[:, :, n] = Y
    X_T[:, :, n] = Y.T

residuals = np.array([])

for i in range(2*3*(nnode-1)):
    r = (X_T[:, :, i] @ (H[:, :, i] @ enlarged_X[:, :, i])) \
    + (g[0,:,i] @ enlarged_X[:,:,i]) \
    + b[0,0,i]
    residuals = np.append(residuals, r)



FTKCL = np.reshape(residuals, (len(residuals), 1))
print(H.shape)
print(g.shape)
print(b.shape)
print(X_T.shape)
print(enlarged_X.shape)
print(nnode)
print(nline)
print(FTKCL.shape)

#--------------- old FTKCL below --------------------

    # ----------- It basically starts here
    #
    # H = np.zeros((2 * 3 * (nnode + nline), 2 * 3* (nnode + nline), 2*3*nline))
    #
    # g = np.zeros((1, 2*3*(nnode+nline), 2*3*nline))
    #
    # b = np.zeros((1, 1, 2*3*nline))
    #
    # for ph in range(0,3):
    #     if ph == 0: #set nominal voltage based on phase
    #         A0 = 1
    #         B0 = 0
    #     elif ph == 1:
    #         A0 = -1/2
    #         B0 = -1 * np.sqrt(3)/2
    #     elif ph == 2:
    #         A0 = -1/2
    #         B0 = np.sqrt(3)/2
    #     for line in range(len(dss.Lines.AllNames())):
    #         for k2 in range(1, len(dss.Circuit.AllBusNames())):
    #             dss.Circuit.SetActiveBus(dss.Circuit.AllBusNames()[k2])
    #             in_lines, out_lines = linelist(dss.Circuit.AllBusNames()[k2]) #get in/out lines of bus
    #             for cplx in range(0,2):
    #                 load_val = d_factor(dss.Circuit.AllBusNames()[k2], cplx)
    #                 gradient_mag = np.array([A0 * ((A0**2+B0**2) ** (-1/2)), B0 * ((A0**2+B0**2) ** (-1/2))]) #some derivatives
    #                 hessian_mag = np.array([[-((A0**2)*(A0**2+B0**2)**(-3/2))+(A0**2+B0**2)**(-1/2), -A0*B0*(A0**2+B0**2)**(-3/2)],
    #                                     [-A0*B0*(A0**2+B0**2)**(-3/2), -((B0**2)*(A0**2+B0**2)**(-3/2))+((A0**2+B0**2)**(-1/2))]])
    #
    #                 #quadratic terms
    #                 H[2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2][2*ph*nline + line*2 + cplx] = -load_val * (beta_Z + (0.5 * beta_I* hessian_mag[0][0]))#check this line for loadval
    #                 H[2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2 + 1][2*ph*nline + line*2 + cplx] = -load_val * (beta_Z + (0.5 * beta_I * hessian_mag[1][1]))
    #                 H[2*(nnode)*ph + 2*k2][2*(nnode)*ph + 2*k2 + 1][2*ph*nline + line*2 + cplx] = -load_val * beta_I * hessian_mag[0][1]
    #                 H[2*(nnode)*ph + 2*k2 + 1][2*(nnode)*ph + 2*k2][2*ph*nline + line*2 + cplx] =  -load_val * beta_I * hessian_mag[0][1]
    #
    #                 for i in range(len(in_lines)): #fill in H for the inlines
    #                     dss.Lines.Name(in_lines[i])
    #                     line_idx = get_line_idx(in_lines[i])
    #                     if cplx == 0: #real residual
    #                         #A_m and C_lm
    #                         H[2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*ph*nline + line*2 + cplx] = 1/2
    #                         H[2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2][2*ph*nline + line*2 + cplx] = 1/2
    #                         #B_m and D_lm
    #                         H[2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*ph*nline + line*2 + cplx] = 1/2
    #                         H[2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1][2*ph*nline + line*2 + cplx] = 1/2
    #                     if cplx == 1: #complex residual
    #                         #A_m, D_lm
    #                         H[2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*ph*nline + line*2 + cplx] = -1/2
    #                         H[2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2][2*ph*nline + line*2 + cplx] = -1/2
    #                         #B_m and C_lm
    #                         H[2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*ph*nline + line*2 + cplx] = 1/2
    #                         H[2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1][2*ph*nline + line*2 + cplx] = 1/2
    #
    #                 for j in range(len(out_lines)): #fill in H for the outlines
    #                     dss.Lines.Name(out_lines[j])
    #                     line_idx = get_line_idx(out_lines[j])
    #                     k = get_bus_idx(dss.Lines.Bus2())
    #                     if cplx == 0:
    #                         #A_m and C_mn
    #                         H[2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*ph*nline + line*2 + cplx] = -1/2
    #                         H[2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2][2*ph*nline + line*2 + cplx] = -1/2
    #                         #B_m and D_mn
    #                         H[2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*ph*nline + line*2 + cplx] = -1/2
    #                         H[2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2 + 1][2*ph*nline + line*2 + cplx] = -1/2
    #                     if cplx == 1:
    #                         #A_m and D_mn
    #                         H[2*(nnode)*ph + 2*k2][2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*ph*nline + line*2 + cplx] = 1/2
    #                         H[2*3*(nnode) + 2*ph*nline + 2*line_idx + 1][2*(nnode)*ph + 2*k2][2*ph*nline + line*2 + cplx] = 1/2
    #                         #C_m and B_mn
    #                         H[2*(nnode)*ph + 2*k2 + 1][2*3*(nnode) + 2*ph*nline + 2*line_idx][2*ph*nline + line*2 + cplx] = -1/2
    #                         H[2*3*(nnode) + 2*ph*nline + 2*line_idx][2*(nnode)*ph + 2*k2 + 1][2*ph*nline + line*2 + cplx] = -1/2
    #
    # #Linear Term
    # for ph in range(0,3):
    #     if ph == 0: #set nominal voltage based on phase
    #         A0 = 1
    #         B0 = 0
    #     elif ph == 1:
    #         A0 = -1/2
    #         B0 = -1 * np.sqrt(3)/2
    #     elif ph == 2:
    #         A0 = -1/2
    #         B0 = np.sqrt(3)/2
    #     for k3 in range(len(dss.Lines.AllNames())):
    #         for k2 in range(len(dss.Circuit.AllBusNames())):
    #             for cplx in range(0,2):
    #                 dss.Lines.Name(dss.Lines.AllNames()[k3])
    #                 bus1 = dss.Lines.Bus1()
    #                 bus2 = dss.Lines.Bus2()
    #                 if bus1 != dss.Circuit.AllBusNames()[k2] or bus2 != dss.Circuit.AllBusNames()[k2]:
    #                     continue
    #
    #                 load_val = d_factor(bus1, cplx)
    #
    #                 #linear terms
    #                 g_temp = np.zeros(len(X))
    #                 g_temp[2*ph*nnode+ 2 * k2] = -load_val* beta_I* ((1/2 * (-2 * A0 * hessian_mag[0][0] - 2 * B0 * hessian_mag[0][1])) \
    #                                        +   gradient_mag[0])
    #                 g_temp[2*ph*nnode+ 2 * k2 + 1] = -load_val * beta_I * ((1/2 * (-2* A0 *hessian_mag[0][1] - 2 * B0 * hessian_mag[1][1])) \
    #                                            +  gradient_mag[1])
    #
    #                 g[0,:,2*nline*ph + 2*k3 + cplx] = g_temp
    #
    # #Constant Term
    # for ph in range(0,3):
    #     if ph == 0: #set nominal voltage based on phase
    #         A0 = 1
    #         B0 = 0
    #     elif ph == 1:
    #         A0 = -1/2
    #         B0 = -1 * np.sqrt(3)/2
    #     elif ph == 2:
    #         A0 = -1/2
    #         B0 = np.sqrt(3)/2
    #     for k3 in range(len(dss.Lines.AllNames())):
    #         for k2 in range(len(dss.Circuit.AllBusNames())):
    #             for cplx in range(0,2):
    #                 dss.Lines.Name(dss.Lines.AllNames()[k3])
    #                 bus1 = dss.Lines.Bus1()
    #                 bus2 = dss.Lines.Bus2()
    #                 if bus1 != dss.Circuit.AllBusNames()[k2] or bus2 != dss.Circuit.AllBusNames()[k2]:
    #                     continue
    #                 load_val = d_factor(bus1, cplx)
    #                 #constant terms
    #                 b_factor = 0
    #                 Sk = dss.CktElement.Powers() #retrieve powers
    #                 if cplx == 1:
    #                     b_factor = dss.Capacitors.kvar() - Sk[1] #depends on if it's real or im residual
    #                 elif cplx == 0:
    #                     b_factor = - Sk[0]
    #
    #                 #i think this is right
    #                 b_temp = -load_val * (beta_S \
    #                     + (beta_I) * (hessian_mag[0][1] * A0 * B0 + (1/2)*hessian_mag[0][0] * A0**2 + (1/2)*hessian_mag[1][1] * B0**2) \
    #                     - beta_I * (A0 * gradient_mag[0] +B0* gradient_mag[1]) \
    #                     + beta_I * (A0**2 + B0**2) ** (1/2)) \
    #                     + b_factor #calculate out the constant term in the residual
    #                 b[0][0][2*nline*ph + 2 *k3+cplx] = b_temp #store the in the b matrix
    #
    #
    # Y = X.reshape(-1, 1) #(72,1)
    # enlarged_X = np.zeros((2*3*(nline+nnode), 1, 2*3*nline))
    # X_T= np.zeros((1, 2*3*(nline+nnode), 2*3*nline))
    # for n in range(2*3*nline):
    #     enlarged_X[:, :, n] = Y
    #     X_T[:, :, n] = Y.T
    #
    # residuals = np.array([])
    #
    # for i in range(2*3*nline):
    #     r = (X_T[:, :, i] @ (H[:, :, i] @ enlarged_X[:, :, i])) \
    #     + (g[0,:,i] @ enlarged_X[:,:,i]) \
    #     + b[0,0,i]
    #     residuals = np.append(residuals, r)
    #
    # FTKCL = np.reshape(residuals, (len(residuals), 1))
