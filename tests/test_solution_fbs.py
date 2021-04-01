# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 21 March 2021
# Test refactor of nr3 module

import numpy as np
import pytest

from circuit_mapper.solution_fbs import SolutionFBS
from circuit_mapper.pretty_print import compare_data_frames
from circuit_mapper.circuit import Circuit

from fbs.fbs.utils import init_from_dss
from fbs.fbs.fbs import fbs, get_solution as get_fbs_solution
import sys
import os
from pathlib import Path

DSS_FILE_DIR = Path('./src/nr3_python/')
OUT_DIR = Path('./tests/test_compare_old_fbs')
OUT_PREFIX = 'FBS_v_FBS_'

GENEROUS = 1e-1
STRICT = 1e-2

def setup_module():
    if not OUT_DIR.exists():
        os.makedirs(OUT_DIR)

@pytest.mark.parametrize(
    "dss_file,tolerance",
    [
        ('IEEE_13_Bus_allwye.dss', GENEROUS),
        ('IEEE_13_Bus_allwye_noxfm_noreg.dss', STRICT),
        ('IEEE_34_Bus_allwye.dss', GENEROUS), # fails
        ('IEEE_34_Bus_allwye_noxfm_noreg.dss', STRICT), # all fail
        ('IEEE_37_Bus_allwye.dss', GENEROUS), # I, sV off
        ('IEEE_37_Bus_allwye_noxfm_noreg.dss', STRICT)
    ],
)
class TestParamtrized:

    def old_fbs_solution(self, dss_file):
        fp = str(Path(DSS_FILE_DIR, dss_file))
        network = init_from_dss(fp)
        return get_fbs_solution(fbs(network))

    def new_fbs_solution(self, dss_file):
        fp = str(Path(DSS_FILE_DIR, dss_file))
        solution = SolutionFBS(fp)
        # solution.maxiter = 1
        solution.solve()
        return solution

    @pytest.mark.parametrize(
        "new_fbs_param, old_fbs_method",
        [
            ('V', 'V_df'),
            ('I', 'I_df'),
            ('sV', 'sV_df'),
        ]
    )
    def test_fbs_v_fbs(self, dss_file, tolerance, new_fbs_param,
                    old_fbs_method):
        """ Test all Files X (V, I, sV). Save log files to OUT_DIR"""
        fp = Path(DSS_FILE_DIR, dss_file)
        Circuit.set_zip_values(np.asarray([0.10, 0.05, 0.85, 0.10, 0.05, 0.85, 0.80]))
        new_fbs_solution = self.new_fbs_solution(dss_file)
        old_fbs_solution = self.old_fbs_solution(dss_file)
        new = new_fbs_solution.get_data_frame(new_fbs_param)
        old = getattr(old_fbs_solution, old_fbs_method)()
        out_file = Path(OUT_DIR, OUT_PREFIX + str(fp.stem) + '_' + new_fbs_param).with_suffix('.out.txt')
        sys.stdout = open(out_file, 'w')
        test = ((new - old).abs().max() <= tolerance).all()
        if test:
            print("TEST PASSED")
        else:
            print("TEST FAILED")
        compare_data_frames(new, old, 'new_fbs', 'old_fbs', new_fbs_param)
        # assert np.testing.assert_allclose(new, old, atol=tolerance)
        assert test

    def test_mapping_FZpu(self, dss_file, tolerance):
        old_network = self.old_fbs_solution(dss_file).network
        new_circuit = self.new_fbs_solution(dss_file).circuit
        # look at FZpu, only for the Lines (not transformers or regs)
        num_pure_lines = new_circuit.lines.num_elements
        FZnew = np.asarray([line.FZpu for line in new_circuit.lines.get_elements()])
        FZold = np.asarray([line.FZpu for line in old_network.get_lines()])
        np.testing.assert_array_equal(FZnew[0: num_pure_lines], FZold[0: num_pure_lines])

    def test_mapping_spu_loads(self, dss_file, tolerance):
        old_network = self.old_fbs_solution(dss_file).network
        new_circuit = self.new_fbs_solution(dss_file).circuit
        spu_new = np.asarray(new_circuit.loads._get_attr_by_idx('spu', 'rows'), dtype=complex)
        # spu_scalar = np.asarray([load.spu for load in old_network.get_loads()])
        spu_old = np.zeros((len(old_network.get_loads()), 3), dtype=complex)
        for load_idx, load in enumerate(old_network.get_loads()):
            for ph_idx, ph in enumerate(load.phases):
                if ph:
                    spu_old[load_idx, ph_idx] = load.spu
        np.testing.assert_equal(spu_new, spu_old)
    
