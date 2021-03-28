# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Test refactor of nr3 and fbs modules

import numpy as np
import opendssdirect as dss
import pytest

from circuit_mapper.circuit import Circuit
from circuit_mapper.solution_nr3 import SolutionNR3
from circuit_mapper.solution import Solution

# current nr3 dependencies
from nr3_python.lib.DSS_parameters import relevant_openDSS_parameters
from nr3_python.lib.basematrices import basematrices
from nr3_python.lib.change_KCL_matrices import change_KCL_matrices
from nr3_python.lib.helper import transformer_regulator_parameters, \
        nominal_load_values, cap_arr, linelist

LOCAL_DIR = 'src/nr3_python/'
DSS_FILE = LOCAL_DIR + 'IEEE_13_Bus_allwye.dss'
# DSS_FILE = LOCAL_DIR + 'IEEE_13_Bus_allwye_noxfm_noreg.dss'

SLACKIDX = 0
VSLACK = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])


@pytest.fixture
def nr3_solution():
    nr3 = SolutionNR3(DSS_FILE)
    nr3.circuit.set_zip_values([1, 0, 0, 1, 0, 0, .8])
    return nr3


@pytest.fixture
def circuit(nr3_solution):
    return nr3_solution.circuit


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
    np.testing.assert_equal(TXnum, circuit.get_tx_idx_matrix())


def test_nr3_DSS_parameters_RXnum(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    np.testing.assert_equal(RXnum, circuit.get_rx_idx_matrix())


def test_nr3_DSS_parameters_PH(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    assert (PH == circuit.buses.get_phase_matrix('cols')).all()


def test_nr3_DSS_parameters_spu(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    assert (spu == circuit.get_spu_matrix()).all()


def test_nr3_DSS_parameters_Z(circuit, nr3_DSS_parameters):
    TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu = \
        nr3_DSS_parameters
    np.testing.assert_equal(aPQ, circuit.get_aPQ_matrix())
    np.testing.assert_equal(aI, circuit.get_aI_matrix())
    np.testing.assert_equal(aZ, circuit.get_aZ_matrix())


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


def test_nr3_linelist(circuit):
    for bus_name in circuit.buses.all_names():
        in_idx = circuit.lines.get_line_list(bus_name, 'in')
        in_lines = linelist(bus_name)[0]
        np.testing.assert_equal(in_idx, in_lines)

    for bus_name in circuit.buses.all_names():
        out_idx = circuit.lines.get_line_list(bus_name, 'out')
        out_lines = linelist(bus_name)[1]
        np.testing.assert_equal(out_idx, out_lines)


def test_nr3_transformer_regulator_params_gain(circuit, xfm_vr_parameters):
    tf_bus, vr_bus, tf_lines, vr_lines, tf_count, vr_no, gain = \
        xfm_vr_parameters
    assert (gain == circuit.voltage_regulators.get_gain_matrix()).all()


def test_nr3_basematrices_XNR(nr3_solution, nr3_basematrices):
    XNR, g_SB, b_SB, G_KVL, b_KVL, H, g, b, H_reg, G_reg = nr3_basematrices
    print_compare("XNR", XNR, nr3_solution.XNR)
    np.testing.assert_equal(XNR, nr3_solution.XNR)


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
    print_compare('H', H, nr3_solution.H)
    np.testing.assert_equal(H, nr3_solution.H)
    assert (g == nr3_solution.g).all()
    assert (b == nr3_solution.b).all()


def print_compare(title, old, new):
    print(title, '-'*70)
    print("old:")
    print(old)
    print("new:")
    print(new)


# NR3 CHANGE KCL TESTS---------------------------------------------------------
def test_change_KCL(nr3_solution):
    H, g, b = nr3_solution.H, nr3_solution.g, nr3_solution.b
    t, der, capacitance = 1, -1, 0
    wpu = nr3_solution.circuit.get_wpu_matrix()
    ckt_H, ckt_b = nr3_solution.change_KCL_matrices(H, g, b, t, der, capacitance, wpu)
    nr3_H, nr3_b = change_KCL_matrices(H, g, b, t, der, capacitance, wpu)
    np.testing.assert_equal(ckt_H, nr3_H)
    np.testing.assert_equal(ckt_b, nr3_b)
