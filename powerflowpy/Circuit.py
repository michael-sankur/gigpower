# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 19 February 2021
# Create Circuit class to mirror a dss Circuit object
# used by Solution objects to solve powerflow

from collections import OrderedDict
import numpy as np
import pandas as pd

import Bus
import Line
from utils import parse_phases


class Circuit():

    def __init__(self, dss, Sbase=10**6):
        """ initialize Circuit from an opendss object's current state"""
        self.Sbase = Sbase
        self._init_buses()
        self._init_lines()

    def _init_buses(self, dss):
        self._all_buses = OrderedDict()  # { bus_name: Bus }
        self._adj_matrix = OrderedDict()  # bus_name: [list of downstream buses]

        all_bus_names = dss.Circuit.AllBusNames()
        self._bus_idx = {b: i for i, b in enumerate(all_bus_names)}

        for node_name in dss.Circuit.AllNodeNames():
            name, phase = node_name.split('.')
            if name not in self._all_buses.keys():
                bus = Bus(name, dss)
                dss.Circuit.SetActiveBus(name)
                bus.Vbase = dss.Bus.kVBase() * 1000
                bus.Ibase = self.Sbase/bus.Vbase
                bus.Zbase = bus.Vbase/bus.Ibase
                self._all_buses[name] = bus
                self._adj_matrix[name] = []
            bus = self._all_buses[name]
            bus.phases.append(phase)

        # iterate through buses to parse phase lists and set bus phase matrix
        bus_ph_matrix = np.zeros((3, len(all_bus_names)), dtype=bool)
        for bus, idx in self._bus_idx.items():
            bus.phases = parse_phases(bus.phases)
            bus_ph_matrix[idx] = bus.phases

        self.bus_ph_df = pd.DataFrame(data=bus_ph_matrix,
                                      index=all_bus_names, columns=['A', 'B', 'C'])

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

        self.ntransformer = 0
        # for transformer_name in transformer_names:
        #     transformer_data = all_transformer_data[transformer_name]
        #     dss.Transformers.Name(transformer_name)

        #     bus1, bus2 = (b.split('.')[0] for b in dss.CktElement.BusNames())
        #     transformer = Transformer(network, (bus1, bus2), transformer_name, transformer_data['NumWindings'])

        #     dss.Transformers.Wdg(1)  # upstream side
        #     if dss.CktElement.NumPhases() == 1:  # 1 phase -> LN voltage,  2 or 3 -> LL voltage
        #         upstream = dss.Transformers.kV()
        #     else:
        #         upstream = dss.Transformers.kV() / (3 ** 0.5)

        #     dss.Transformers.Wdg(2)  # downstream side
        #     if dss.CktElement.NumPhases() == 1:  # 1 phase -> LN voltage,  2 or 3 -> LL voltage
        #         downstream = dss.Transformers.kV()
        #     else:
        #         downstream = dss.Transformers.kV() / (3 ** 0.5)
        #     transformer.rated_voltages = (upstream, downstream)

        #     transformer.phases = network.nodes[bus1].phases
        #     transformer.conn = 'delta' if transformer_data['IsDelta'] else 'wye'

        #     transformer.kV = transformer_data['kV']
        #     transformer.kVA = transformer_data['kVA']




    #NODE indexing
    TXnum = np.zeros((nline), dtype='int')
    RXnum = np.zeros((nline), dtype='int')
    TXnode = [None]*nline
    RXnode = [None]*nline  # name of outgoing bus on line


    #TXnum and RXnum
    for line in range(len(dss.Lines.AllNames())):
        dss.Lines.Name(dss.Lines.AllNames()[line])  # set the line
        bus1 = dss.Lines.Bus1()
        bus2 = dss.Lines.Bus2()
        pattern = r"(\w+)."  # this appears to wrok

        TXnode[line] = re.findall(pattern, bus1)[0]
        RXnode[line] = re.findall(pattern, bus2)[0]
        TXnum[line] = dss.Circuit.AllBusNames().index(TXnode[line])
        RXnum[line] = dss.Circuit.AllBusNames().index(RXnode[line])

    #TF
    for line in range(len(line_idx_tf)):
        print('is there a trnasformer?')
        lineidx = line_idx_tf[line]
        TXnode[lineidx] = dss.Circuit.AllBusNames()[tf_bus[0, line]]  # bus name
        RXnode[lineidx] = dss.Circuit.AllBusNames()[tf_bus[1, line]]
        TXnum[lineidx] = tf_bus[0, line]  # bus index
        RXnum[lineidx] = tf_bus[1, line]
    #VR in
    for line in range(len(line_in_idx_vr)):
        print('is there a voltage regulator?')
        lineinidx = line_in_idx_vr[line]

        TXnode[lineinidx] = dss.Circuit.AllBusNames()[vr_bus[0, line]]
        RXnode[lineinidx] = dss.Circuit.AllBusNames()[vr_bus[1, line]]
        TXnum[lineinidx] = vr_bus[0, line]
        RXnum[lineinidx] = vr_bus[1, line]
    #VR out
    for line in range(len(line_out_idx_vr)):
        print('is there a voltage regulator part 2?')
        lineoutidx = line_out_idx_vr[line]
        TXnode[lineoutidx] = dss.Circuit.AllBusNames()[vr_bus[0, line]]
        RXnode[lineoutidx] = dss.Circuit.AllBusNames()[vr_bus[1, line]]
        TXnum[lineoutidx] = vr_bus[0, line]
        RXnum[lineoutidx] = vr_bus[1, line]

    #spu, apq, ai, az
    spu = np.zeros((3, nnode))
    ppu = np.zeros((3, nnode))
    qpu = np.zeros((3, nnode))
    aPQ = np.zeros((3, nnode))
    aI = np.zeros((3, nnode))
    aZ = np.zeros((3, nnode))

    if t == -1:
        var = 1
    else:
        var = (1 + 0.1*np.sin(2*np.pi*0.01*t))

    for n in range(len(dss.Loads.AllNames())):  # go through the loads
        dss.Loads.Name(dss.Loads.AllNames()[n])  # set the load
        load_phases = [0, 0, 0]  # instantiate load phases as all non-existent
        pattern = r"(\w+)\."
        load_bus = re.findall(pattern, dss.CktElement.BusNames()[0])  # determine bus name
        knode = dss.Circuit.AllBusNames().index(load_bus[0])  # match busname to index
        for ph in range(0, 3):
            pattern = r"\.%s" % (str(ph + 1))
            m = re.findall(pattern, dss.CktElement.BusNames()[0])
            if m:  # if phase exists for load
                load_phases[ph - 1] = 1
                aPQ[ph, knode] = 1
                aZ[ph, knode] = 0
                ppu[ph, knode] = dss.Loads.kW() * 1e3 * var / Sbase  # check these lines later
                qpu[ph, knode] = dss.Loads.kvar() * 1e3 * var / Sbase
        if sum(load_phases) > 1:  # if it is a multiphase load
            for ph in range(0, 3):
                ppu[ph, knode] /= sum(load_phases)
                qpu[ph, knode] /= sum(load_phases)
    spu = (ppu + 1j * qpu)

    #cappu, wpu, vvcpu
    cappu = np.zeros((3, nnode))

    def cap_dict():
        for n in range(len(dss.Capacitors.AllNames())):
            dss.Capacitors.Name(dss.Capacitors.AllNames()[n])
            cap_data = dss.CktElement.BusNames()[0].split('.')

            idxbs = dss.Circuit.AllBusNames().index(cap_data[0])
            for ph in range(1, len(cap_data)):
                cappu[int(cap_data[ph]) - 1, idxbs] += dss.Capacitors.kvar() * 1e3 / Sbase / (len(cap_data) - 1)

    wpu = np.zeros((3, nnode))
    vvcpu = np.zeros((3, nnode))

    return TXnum, RXnum, PH, spu, aPQ, aZ, aI, cappu, wpu, vvcpu


