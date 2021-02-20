# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: October 2020
# Make sure that fbs module fails on a network that is not radial.

from fbs.utils import init_from_dss
from fbs.fbs import fbs
import pytest

dss_file_cycle = 'fbs/tests/non_radial_networks/06node_contains_cycle.dss'
dss_file_multiple_parents = 'fbs/tests/non_radial_networks/06node_multiple_parents.dss'


@pytest.fixture
def get_cycle_network():
    """setup test network from dss_file"""
    return init_from_dss(dss_file_cycle)


def test_fail_on_cycle(get_cycle_network) -> None:
    """Make sure that topo_sort raises a Value Error on a network that contains a cycle"""
    with pytest.raises(ValueError):
        fbs(get_cycle_network)


def test_fail_on_multiple_parents() -> None:
    """
    Make sure that init_from_dss raises a Value Error on a network
    with a node that has more than one parent
    """
    with pytest.raises(ValueError):
        init_from_dss(dss_file_multiple_parents)
