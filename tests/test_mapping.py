# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Test refactor of nr3 and fbs modules

import numpy as np
import opendssdirect as dss
import pytest

from circuit import Circuit

# previous nr3 mapping functions
from nr3.lib.compute_vecmat import compute_vecmat
from nr3.lib.relevant_openDSS_parameters import relevant_openDSS_parameters

LOCAL_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
DSS_FILE = LOCAL_DIR + 'IEEE_13_Node/IEEE_13_Bus_allwye_noxfm_noreg.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])


# @pytest.fixture
def circuit():
    """ map Circuit object once for use in all tests """
    dss.run_command('Redirect ' + DSS_FILE)
    return Circuit(dss)


# NR3 TESTS---------------------------------------------------------------------
def test_nr3_relevant_open_DSS_parameters():
    dss.run_command('Redirect ' + DSS_FILE)
    circuit = Circuit(dss)

    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        relevant_openDSS_parameters(DSS_FILE, -1)
    assert (TXnum == circuit.get_lines_tx_idx_matrix()).all()
    assert (RXnum == circuit.get_lines_rx_idx_matrix()).all()
    # assert PH == circuit.buses.get_phase_matrix()
    # assert spu == circuit.buses.get_spu_matrix()
    # assert aPQ == circuit.buses.get_aPQ_matrix()
    # assert aI == circuit.buses.get_aZ_matrix()
    # assert aZ == circuit.buses.get_aI_matrix()
    # assert cappu == circuit.buses.get_cappu_matrix()
    # assert wpu == circuit.buses.get_wpu_matrix()
    # assert vvcpu == circuit.buses.get_vvcpu_matrix()


# def test_nr3_compute_vecmat(circuit):
#     X, g_SB, b_SB, G_KVL, b_kvl, H_reg, G_reg = compute_vecmat(circuit.get_XNR_mat(), DSS_FILE, VSLACK)
#     assert X == circuit.calc_nr3_X()
#     assert g_SB == circuit.calc_nr3_g_SB()
#     assert b_SB == circuit.calc_nr3_b_SB()
#     assert G_KVL == circuit.calc_nr3_G_KVL()
#     assert b_KVL == circuit.calc_nr3_b_kvl()
#     assert H_reg == circuit.calc_nr3_H_reg()
#     assert G_reg == circuit.calc_nr3_G_reg()
