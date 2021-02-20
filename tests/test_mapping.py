# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Test refactor of nr3 and fbs modules

import numpy as np
import opendssdirect as dss
import pytest

# previous nr3 mapping functions
from nr3.lib.compute_vecmat import compute_vecmat
from nr3.lib.relevant_openDSS_parameters import relevant_openDSS_parameters


# current mapper, the Circuit object
import Circuit

LOCAL_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
DSS_FILE = LOCAL_DIR + 'IEEE_13_Node/IEEE_13_Bus_allwye_noxfm_noreg.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])


@pytest.fixture
def circuit():
    """ map Circuit object once for use in all tests """
    dss.run_command('Redirect ' + DSS_FILE)
    return Circuit(dss)


# NR3 TESTS---------------------------------------------------------------------
def test_nr3_relevant_open_DSS_parameters(circuit):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = relevant_openDSS_parameters(DSS_FILE)
    assert TXnum == circuit.get_line_TXnum()
    assert RXnum == circuit.get_line_RXnum()
    assert PH == circuit.get_bus_ph_mat()
    assert spu == circuit.get_bus_spu_mat()
    assert aPQ == circuit.get_bus_aPQ_mat()
    assert aI == circuit.get_bus_aZ_mat()
    assert aZ == circuit.get_bus_aI_mat()
    assert cappu == circuit.get_bus_cappu_mat()
    assert wpu == circuit.get_bus_wpu_mat()
    assert vvcpu == circuit.get_bus_vvcpu_mat()


def test_nr3_compute_vecmat(circuit):
    X, g_SB, b_SB, G_KVL, b_kvl, H_reg, G_reg = compute_vecmat(circuit.get_XNR_mat(), DSS_FILE, VSLACK)
    assert X == circuit.calc_nr3_X()
    assert g_SB == circuit.calc_nr3_g_SB()
    assert b_SB == circuit.calc_nr3_b_SB()
    assert G_KVL == circuit.calc_nr3_G_KVL()
    assert b_KVL == circuit.calc_nr3_b_kvl()
    assert H_reg == circuit.calc_nr3_H_reg()
    assert G_reg == circuit.calc_nr3_G_reg()
