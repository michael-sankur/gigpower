# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Circuit class to mirror a dss Circuit object
# used by Solution objects to solve powerflow

from collections import OrderedDict
import numpy as np
import pandas as pd

import BusGroup, CapacitorGroup, LineGroup, LoadGroup, TransformerGroup, VoltageRegulatorGroup



class Circuit():

    def __init__(self, dss, Sbase=10**6):
        """ initialize Circuit from an opendss object's current state"""
        self.Sbase = Sbase
        self._init_buses()
        self._init_lines()

    def _init_buses(self, dss):
        self._all_buses = OrderedDict()  # { bus_name: Bus }
        self._adj_matrix = OrderedDict()  # bus_name: [list of downstream buses]



    def _init_lines(self, dss):
        self.nline = 0
        all_lines_data = dss.utils.lines_to_dataframe().transpose()
        line_codes = all_lines_data.keys()
        # TODO: set line base values based on node Vbase
        for line_code in line_codes:
            line_data = all_lines_data[line_code]
            # Handle bus names that include phases
            if "." in line_data['Bus1']:
                tx, *tx_phases = line_data['Bus1'].split('.')
                rx, *rx_phases = line_data['Bus2'].split('.')
                if tx_phases != rx_phases:
                    raise ValueError(f'Tx phases do not match Rx phases for line {line_code}')
            else:  # Otherwise, define the line for all 3 phases.
                tx, rx = line_data['Bus1'], line_data['Bus2']
                tx_phases = ['1', '2', '3']

            line = Line(line_code, tx, rx)
            line.phases = parse_phases(tx_phases)
            tx_node = network.nodes.get(tx)

            # parse line attributes from dss line data
            line.length = line_data['Length']
            # TODO: HERE 2/19

            fz_mult = 1 / tx_node.Zbase * line.length
            line.FZpu = get_Z(line_data, line.phases, fz_mult)
            network.lines[self.key] = self
            tx, rx = self.key
            if rx not in network.adj[tx]:
                network.adj[tx].append(rx)  # add this as a Line in the adjacency list
            # store a pointer to tx on rx.parent
            if network.nodes[rx].parent:
                raise ValueError(f"Error when processing line {self.name}. Node {rx} already has a parent.")
            network.nodes[rx].parent = network.nodes[tx]

        pass


    def _init_transformers(self, dss):
        """ map Transformer objects. Save to network in an Ordered Dict"""
        all_transformer_data = dss.utils.transformers_to_dataframe().transpose()
        transformers = dss.Transformers.AllNames()
        # exclude the transformers that are paired with voltage regulators
        vr_transformers = [r.transformer_name for r in network.voltageRegulators.values()]
        transformer_names = [n for n in transformers if n not in vr_transformers]




