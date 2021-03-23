# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Test refactor of nr3 and fbs modules

import numpy as np
import opendssdirect as dss
import pytest

from circuit import Circuit
from solution_nr3 import SolutionNR3
from solution import Solution

# current nr3 dependencies
import sys
#sys.path.append('/Users/elainelaguerta/Dropbox/LBNL/LinDist3Flow/20180601/PYTHON')
sys.path.append("C:\\Users\\kathl\\Desktop\\GitHub\\LinDist3Flow\\20180601\\PYTHON")

from lib.DSS_parameters import relevant_openDSS_parameters
from lib.basematrices import basematrices
from lib.helper import transformer_regulator_parameters, nominal_load_values, cap_arr

LOCAL_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
LOCAL_DIR = "C:\\Users\\kathl\\Desktop\\GitHub\\LinDist3Flow\\20180601\\PYTHON"
DSS_FILE = LOCAL_DIR + '\\IEEE_13_Bus_allwye.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])


@pytest.fixture
def circuit():
    """ map Circuit object once for use in all tests """
    dss.run_command('Redirect ' + DSS_FILE)
    return Circuit(dss)


@pytest.fixture
def nr3_solution():
    return SolutionNR3(DSS_FILE)


@pytest.fixture
def nr3_DSS_parameters():
    return relevant_openDSS_parameters(DSS_FILE, [])


@pytest.fixture
def xfm_vr_parameters():
    return transformer_regulator_parameters()


@pytest.fixture
def nr3_basematrices():
    slack_idx = 0
    Vslack = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])
    V0, I0 = None, None
    return basematrices(slack_idx, Vslack, V0, I0)


# NR3 TESTS---------------------------------------------------------------------
def test_nr3_DSS_parameters_TXnum(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    print_compare("TX", TXnum, circuit.get_tx_idx_matrix())
    assert (TXnum == circuit.get_tx_idx_matrix()).all()


def test_nr3_DSS_parameters_RXnum(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    print_compare("RX", RXnum, circuit.get_rx_idx_matrix())
    assert (RXnum == circuit.get_rx_idx_matrix()).all()


def test_nr3_DSS_parameters_PH(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    print_compare("PH", PH, circuit.buses.get_phase_matrix())
    assert (PH == circuit.buses.get_phase_matrix()).all()


def test_nr3_DSS_parameters_spu(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    print_compare("spu", spu, circuit.get_spu_matrix())
    assert (spu == circuit.get_spu_matrix()).all()


def test_nr3_DSS_parameters_Z(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    assert (aPQ == circuit.get_aPQ_matrix()).all()
    assert (aI == circuit.get_aI_matrix()).all()
    assert (aZ == circuit.get_aZ_matrix()).all()


def test_nr3_DSS_parameters_cappu(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    print_compare("cappu", cappu, circuit.get_cappu_matrix())
    assert (cappu == circuit.get_cappu_matrix()).all()


def test_nr3_DSS_parameters_wpu_vvcpu(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    assert (wpu == circuit.get_wpu_matrix()).all()
    assert (vvcpu == circuit.get_vvcpu_matrix()).all()


def test_nr3_transformer_regulator_params_tf_bus(circuit, xfm_vr_parameters):
    tf_bus, vr_bus, tf_lines, vr_lines, tf_count, vr_no, gain = \
        xfm_vr_parameters
    print_compare('tf_bus', tf_bus, circuit.transformers.get_bus_ph_matrix())
    assert (tf_bus == circuit.transformers.get_bus_ph_matrix()).all()


def test_nr3_transformer_regulator_params_vr_bus(circuit, xfm_vr_parameters):
    tf_bus, vr_bus, tf_lines, vr_lines, tf_count, vr_no, gain = \
        xfm_vr_parameters
    print_compare('vr_bus', vr_bus, circuit.voltage_regulators.get_bus_ph_matrix())
    assert (vr_bus == circuit.voltage_regulators.get_bus_ph_matrix()).all()


def test_nr3_transformer_regulator_params_counts(circuit, xfm_vr_parameters):
    tf_bus, vr_bus, tf_lines, vr_lines, tf_count, vr_no, gain = \
        xfm_vr_parameters
    assert tf_lines == circuit.transformers.get_num_lines_x_ph()
    assert vr_lines == circuit.voltage_regulators.get_num_lines_x_ph()
    assert tf_count == circuit.transformers.num_elements
    assert vr_no == circuit.voltage_regulators.num_elements


def test_nr3_transformer_regulator_params_gain(circuit, xfm_vr_parameters):
    tf_bus, vr_bus, tf_lines, vr_lines, tf_count, vr_no, gain = \
        xfm_vr_parameters
    assert (gain == circuit.voltage_regulators.get_gain_matrix()).all()


def test_nr3_basematrices_XNR(nr3_solution, nr3_basematrices):
    XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg = nr3_basematrices
    print_compare("XNR", XNR, nr3_solution.XNR)
    assert (XNR == nr3_solution.XNR).all()


def test_nr3_basematrices_SB(nr3_solution, nr3_basematrices):
    XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg = nr3_basematrices
    assert (g_SB == nr3_solution.g_SB).all()
    assert (b_SB == nr3_solution.b_SB).all()


def test_nr3_basematrices_KVL(nr3_solution, nr3_basematrices):
    tolerance = 1e-6
    XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg = nr3_basematrices
    assert (b_KVL == nr3_solution.b_KVL).all()
    assert (abs(G_KVL - nr3_solution.G_KVL) <= tolerance).all()


def test_nr3_basematrices_KVL_regs(nr3_solution, nr3_basematrices):
    tolerance = 1e-6
    XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg = nr3_basematrices
    assert (H_reg == nr3_solution.H_reg).all()
    assert (G_reg == nr3_solution.G_reg).all()


def test_nr3_helpers(circuit):
    dsskw, dsskvar = nominal_load_values(-1)
    cappu = cap_arr()
    assert (circuit.loads.get_ppu_matrix == dsskw).all
    assert (circuit.loads.get_qpu_matrix == dsskvar).all
    assert (circuit.capacitors.get_cappu_matrix == cappu).all


def test_nr3_basematrices_H_g_b(nr3_solution, nr3_basematrices):
    XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg = nr3_basematrices
    assert (H == nr3_solution.H).all()
    assert (g == nr3_solution.g).all()
    assert (b == nr3_solution.b).all()


def print_compare(title, old, new):
    print(title, '-'*70)
    print("old:")
    print(old)
    print("new:")
    print(new)

# NR3 CHANGE KCL TESTS---------------------------------------------------------------------

def test_change_KCL(circuit, H, g, b, t, der, capacitance, wpu):
    ckt_H, ckt_b = circuit.change_KCL_matrices( H, g, b, t, der, capacitance, wpu)
    nr3_H, nr3_b = change_KCL_matrices(H, g, b, t, der, capacitance, wpu)
    
    assert(ckt_H == nr3_H)
    assert(ckt_b == nr3_b)