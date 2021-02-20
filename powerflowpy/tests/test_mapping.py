# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Test refactor of nr3 and fbs modules

import numpy as np
import opendssdirect as dss
import pytest

# previous nr3 mapping functions
from nr3.lib.compute_NR3FT_vectorized import compute_NR3FT_vectorized as ft1
from nr3.lib.compute_NR3JT_vectorized import compute_NR3JT_vectorized as jt1
from nr3.lib.compute_vecmat import compute_vecmat
from nr3.lib.relevant_openDSS_parameters import relevant_openDSS_parameters


# import current mapping functions
import Circuit

LOCAL_DIR = '/Users/elainelaguerta/Dropbox/LBNL/python-powerflow/'
DSS_FILE = LOCAL_DIR + 'IEEE_13_Node/IEEE_13_Bus_allwye_noxfm_noreg.dss'


@pytest.fixture
def get_dss_instance():
    """set-up dss instance for dss_file"""
    dss.run_command('Redirect ' + DSS_FILE)
    return None


def test_prev_nr3_mapping():
    """ make sure new mapping is backwards compatible with nr3 functions"""
    pass
