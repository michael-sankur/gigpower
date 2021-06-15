import numpy as np


def map_output(circuit, XNR):

    nline = circuit.lines.num_elements
    nnode = circuit.buses.num_elements

    TXnum = circuit.get_tx_idx_matrix()
    RXnum = circuit.get_rx_idx_matrix()
    PH = circuit.buses.get_phase_matrix('cols')
    spu = circuit.get_spu_matrix()
    APQ = circuit.get_aPQ_matrix()
    AZ = circuit.get_aZ_matrix()
    AI = circuit.get_aI_matrix()
    cappu = circuit.get_cappu_matrix()
    wpu = circuit.get_wpu_matrix()

    # remap XNR to VNR, INR, STXNR, SRXNR, iNR, sNR
    VNR = np.zeros((3, nnode), dtype='complex')
    for ph in range(3):
        for k1 in range(nnode):
            VNR[ph, k1] = XNR[2*ph*nnode + 2*k1] + 1j*XNR[2*ph*nnode + 2*k1+1]
            if np.abs(VNR[ph, k1].real) <= 1e-12:
                VNR[ph, k1] = 0 + VNR[ph, k1].imag
            if np.abs(VNR[ph,k1].imag) <= 1e-12:
                VNR[ph, k1] = VNR[ph, k1].real + 0
    VNR[PH == 0] = 0

    XNR = XNR[2*3*nnode:]
    INR = np.zeros((3,nline), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nline):
            INR[ph,k1] = XNR[2*ph*nline + 2*k1] + 1j*XNR[2*ph*nline + 2*k1+1]
            if np.abs(INR[ph,k1].real) <= 1e-12:
                INR[ph,k1] = 0 + INR[ph,k1].imag
            if np.abs(INR[ph,k1].imag) <= 1e-12:
                INR[ph,k1] = INR[ph,k1].real + 0

    # STXNR_n^phi = V_m^phi (I_mn^phi)^*
    # SRXNR_n^phi = V_n^phi (I_mn^phi)^*
    STXNR = np.zeros((3,nline), dtype='complex')
    SRXNR = np.zeros((3,nline), dtype='complex')
    for ph in range(0,3):
        for k1 in range(0,nline):
            STXNR[ph,k1] = VNR[ph,TXnum[k1]]*np.conj(INR[ph,k1])
            if np.abs(STXNR[ph,k1].real) <= 1e-12:
                STXNR[ph,k1] = 0 + STXNR[ph,k1].imag
            if np.abs(STXNR[ph,k1].imag) <= 1e-12:
                STXNR[ph,k1] = STXNR[ph,k1].real + 0
            SRXNR[ph,k1] = VNR[ph,RXnum[k1]]*np.conj(INR[ph,k1])
            if np.abs(SRXNR[ph,k1].real) <= 1e-12:
                SRXNR[ph,k1] = 0 + SRXNR[ph,k1].imag
            if np.abs(SRXNR[ph,k1].imag) <= 1e-12:
                SRXNR[ph,k1] = SRXNR[ph,k1].real + 0

    sNR = np.zeros((3,nnode), dtype='complex')
    iNR = np.zeros((3,nnode), dtype='complex')
    # Totdal node loads
    sNR = spu*(APQ + AI*np.abs(VNR) + AZ*np.abs(VNR)**2) - 1j*cappu.real + 1j*wpu
    sNR[PH == 0] = 0
    for ph in range(0,3):
        for k1 in range(0,nnode):
            if np.abs(sNR[ph,k1].real) <= 1e-12:
                sNR[ph,k1] = 0 + sNR[ph,k1].imag
            if np.abs(sNR[ph,k1].imag) <= 1e-12:
                sNR[ph,k1] = sNR[ph,k1].real + 0

    # Total node current
    iNR[PH != 0] = np.conj(sNR[PH != 0]/VNR[PH != 0])
    iNR[PH == 0] = 0
    for ph in range(0,3):
        for k1 in range(0,nnode):
            if np.abs(iNR[ph,k1].real) <= 1e-12:
                iNR[ph,k1] = 0 + iNR[ph,k1].imag
            if np.abs(iNR[ph,k1].imag) <= 1e-12:
                iNR[ph,k1] = iNR[ph,k1].real + 0

    return VNR, INR, STXNR, SRXNR, iNR, sNR