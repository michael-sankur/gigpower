# Smoke tests for pre-requisites to running fbs.
# These are: reading network attributes from opendss, and topo-sort

from powerflowpy.utils import init_from_dss
from powerflowpy.fbs import topo_sort

dss_file = 'powerflowpy/tests/test_cases_dss/compare_opendss_05node_threephase_unbalanced_oscillation_03.dss'

def test_init_from_dss() -> None:
    """ Compare opendss.circuit to our Network object"""
    pass
    # number of nodes
    # number of lines

    # nodes x phases

    # lines x phases

def test_topo_sort() -> None:
    network = init_from_dss(dss_file)
    print(topo_sort(network))
    assert True
    # topo sort fails on a network that has a cycle
    # every node precedes its children in the ordering

def test_network_topo() -> None:
    # check that the first node of topo sort is the source node in the dss file
    # check that opendss' assesment of cycles is the same as toposort
    pass
