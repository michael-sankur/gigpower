import numpy as np
import opendssdirect as dss
import sys
from typing import Iterable, List
from network import Node, Line, Load, Controller, Capacitor

def init_from_dss(network: Network, dss_fp: str) -> None:
    """define Network attributes from a dss file"""
    dss.run_command('Redirect ' + dss_fp)
    # set base values
    #TODO: figure out how to get Vbase, units and Sbase, units from opendss

    # make Nodes
    for node_name in dss.Circuit.AllNodeNames():
        name, phase  = node_name.split('.')
        # get the node corresponding to this name, or make a new one
        if name not in network.nodes.keys():
            network.nodes[name] = Node(name)
        node = network.nodes[name]
        phase_idx = int(phase) - 1 # shift to 0-indexing
        self.phases[phase_idx] = 1 # add this phase to the node
    
    #make Lines
    line_codes = dss.LineCodes.AllNames()
    lengths = [i() for i in dss.utils.Iterator(dss.Lines, 'Length')]
    line_zip = zip(line_codes, lengths)
    
    for line_code,length in line_zip:
        tx, rx = line_code.split('_')
        if (tx, rx) not in network.lines.keys():
            network.lines[(tx,rx)] = Line((tx,rx))
        line = network.lines[(tx,rx)]
        line.length = length
    # TODO: figure out how to get 3x3 impedance per unit matrix. 
    # probably something using dss.Lines, 'XMatrix' and 'RMatrix'
    # TODO: handle unit conversions

    # make Loads
    load_names = dss.Loads.AllNames()
    for load_name in load_names:
        # TODO: handle multiple loads on a node 
        node_name, phase_char, load_idx = load_name.split('_')[1:]
        try:
            node = network.nodes[node_name]
        except KeyError ("Node assigned to load not yet defined.")
        load = Load(load_name + '_' + load_idx)
        load.phases[get_phase_idx(phase_char)] = 1 # indicate this load's phase
        node.load = load #assign this load to its node
    #TODO: get aPQ, aI, ppu, qpu, spu for each load
    #TODO: figure out how to get connection information (wye or delta). 
    # Is this info determined by Capacitors?
    # Or can we use Loads.IsDelta?

    # make Controllers
    # TODO: implement this. No idea how opendssdirect maps this info. Which class is it even?

    # make Capacitors
    cap_names = dss.Capacitors.AllNames()
    for cap_name in 

def get_phase_idx(phase_char: str) -> int:
    """
    helper function to turn a phase letter into an index, where 'a' = 0
    """
    return ord(phase_char.lower()) - ord('a') 