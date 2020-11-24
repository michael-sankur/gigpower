# Smoke tests for pre-requisites to running fbs.

from fbs.utils import init_from_dss
from fbs.fbs import topo_sort
import pytest
import opendssdirect as dss
from fbs.network import Transformer

# dss_file = 'fbs/tests/06n3ph_unbal/06node_threephase_unbalanced.dss'
# dss_file = 'fbs/tests/IEEE_13_bus/IEEE_13_Bus_allwye_noxfm_noreg.dss'
dss_file = 'fbs/tests/IEEE_13_bus/IEEE_13_Bus_allwye.dss'
# dss_file = 'fbs/tests/IEEE_13_bus/IEEE_13_Bus_original.dss'

@pytest.fixture
def get_network():
    """setup test network from dss_file"""
    return init_from_dss(dss_file)

@pytest.fixture
def get_dss_instance():
    """set-up dss instance for dss_file"""
    dss.run_command('Redirect ' + dss_file)
    return None

def test_init_from_dss(get_network) -> None:
    """ Compare opendss.circuit to our Network object"""
    # check that number of dss nodes equals our network's node x phases
    network = get_network
    dss_nodes = len(dss.Circuit.AllNodeNames())
    network_nodes = sum([node.phases.count(True) for node in network.get_nodes()])
    node_check = dss_nodes == network_nodes
    # check that number of lines match
    dss_lines = len(dss.Lines.AllNames())

    # count the number of lines in the network that do not represent transformers
    network_lines = len(network.lines) - len(network.transformers)
    # subtract the synthetic upstream lines for voltage regulators
    vr_lines = 0
    for line in network.get_lines():
        if line.voltageRegulators:
            vr_lines += 1
    network_lines -= vr_lines

    # check that the network adjacency list aligns with the number of lines
    network_edges = sum([len(node_list) for node_list in network.adj.values()])
    # subtract out transformers and upstream lines
    network_edges = network_edges - len(network.transformers) - vr_lines
    line_check = (dss_lines == network_lines) and (dss_lines == network_edges)
    assert node_check and line_check

    # check that the number of loads match
    dss_loads = len(dss.Loads.AllNames())
    network_loads = len(network.loads)
    load_check = dss_loads == network_loads

    # check that list of loads on each node matches opendss
    load_list_check = True
    for node in network.get_nodes():
        dss.Circuit.SetActiveBus(node.name)
        if len(dss.Bus.LoadList()) != len(node.loads):
            load_list_check = False

    # check that the number of capacitors match
    dss_capacitors = len(dss.Capacitors.AllNames())
    network_capacitors = len(network.capacitors)
    capacitor_check = dss_capacitors == network_capacitors

    # check that number of transformers match
    dss_transformers = len(dss.Transformers.AllNames())
    network_transformers = len(network.transformers) + len(network.voltageRegulators)
    transformer_check = dss_transformers == network_transformers
    assert node_check and line_check and load_check and load_list_check
    assert capacitor_check and transformer_check

def test_topo_sort() -> None:
    network = init_from_dss(dss_file)
    # if any node is preceded by one of its children assert False
    sort = topo_sort(network)
    for idx,node in enumerate(sort):
        children = network.adj[node]
        for child in children:
            if sort.index(child) < idx:
                assert False

def test_network_topo() -> None:
    # check that the first node of topo sort is the source node in the dss file
    # check that opendss' assesment of cycles is the same as toposort
    # TODO: write me!
    pass
