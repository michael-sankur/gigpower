% Michael Sankur - msankur@lbl.gov
% 2018.06.01

clc, clear, close all

%% Load all network data

% Read from csv file, input first 4 columns as strings (general case)
% Change file path as needed
fp = ['/home/michael/Dropbox/Unbalanced LinDistflow/20180601/NETWORKS/'];
fn = 'ieee_13node_balance.txt';
fn = '05node_fullphase_radial.txt';
fn = '03node_singlephase_radial_example.txt';
fn = '03node_singlephase_mesh_example.txt';
fn = 'ieee_34node.txt';



linearray = {};
fid = fopen([fp fn]);
tline = fgetl(fid);
while ischar(tline)
    linearray{end+1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

% Read lines and store in separate arrays for base values, nodes, lines,
% configurations, loads, capacitors, cons, vvc
baselinearray = {};
nodelinearray = {};
linelinearray = {};
configlinearray = {};
loadlinearray = {};
capacitorlinearray = {};
controllerlinearray = {};
vvclinearray = {};
for k1 = 1:size(linearray,2)
    templinearray = strsplit(linearray{k1},' ');
    if strcmp('base', templinearray{1})
        baselinearray{end+1} = linearray{k1};
    end
    if strcmp('node', templinearray{1})
        nodelinearray{end+1} = linearray{k1};
    end
    if strcmp('line', templinearray{1})
        linelinearray{end+1} = linearray{k1};
    end
    if strcmp('config', templinearray{1})
        configlinearray{end+1} = linearray{k1};
    end
    if strcmp('load', templinearray{1})
        loadlinearray{end+1} = linearray{k1};
    end
    if strcmp('capacitor', templinearray{1})
        capacitorlinearray{end+1} = linearray{k1};
    end
    if strcmp('controller', templinearray{1})
        controllerlinearray{end+1} = linearray{k1};
    end
    if strcmp('vvc', templinearray{1})
        vvclinearray{end+1} = linearray{k1};
    end
end

%% Base Values

% Base values [V], [VAr], [A], [Ohm]

% Iterate through base value lines
for k1 = 1:size(baselinearray,2)
    temp1 = strsplit(baselinearray{k1},' ');
    temp1 = temp1(2:end);
    
    for k2 = 1:length(temp1)        
        temp2 = strsplit(temp1{k2},'=');
        
        % Voltage base value [V]
        if strcmp(temp2{k2},'Vbase')
            Vbase = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('unit',temp2{1})
                if strcmp(temp2{2},'V')
                    Vbase = Vbase;
                end
                if strcmp(temp2{2},'kV')
                    Vbase = 1e3*Vbase;
                end
                if strcmp(temp2{2},'MW')
                    Vbase = 1e6*Vbase;
                end
            end
        end

        % Power base value [VAr]
        if strcmp(temp2{k2},'Sbase')
            Sbase = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('unit',temp2{1})
                if strcmp(temp2{2},'VAr')
                    Sbase = Sbase;
                end
                if strcmp(temp2{2},'kVAr')
                    Sbase = 1e3*Sbase;
                end
                if strcmp(temp2{2},'MVAr')
                    Sbase = 1e6*Sbase;
                end
            end
        end
    end
    
end

% Calculate current and impedance base values
Ibase = Sbase/Vbase; % [A]
Zbase = Vbase/Ibase; % [Ohm]

base.Vbase = Vbase;
base.Sbase = Sbase;
base.Ibase = Ibase;
base.Zbase = Zbase;

% Place base values in network object
% network.Vbase = Vbase;
% network.Sbase = Sbase;
% network.Ibase = Ibase;
% network.Zbase = Zbase;

%% Nodes

% List of possible phases for nodes and corresponding matrix
phlist = {'a','b','c','ab','bc','ac','abc'};
PHmat = [...
1 0 0 1 0 1 1;
0 1 0 1 1 0 1;
0 0 1 0 1 1 1];

% Iterate through node in configuration file. Nodes are added in order of
% their placement in the configuration file. The slacknode should be the
% first node in the configuration file (for the time being)
for k1 = 1:size(nodelinearray,2)
    temp1 = strsplit(nodelinearray{k1},' ');
    temp1 = temp1(2:end);
    
    for k2 = 1:length(temp1)
        temp2 = strsplit(temp1{k2},'=');
        
        % Node names
        if strcmp(temp2{1},'nodename')
            nodes.nodelist(k1) = temp2(2);
        end

        % Node phase list and node phase matrix
        if strcmp(temp2{1},'phases')
            nodes.phases{k1} = temp2{2};
            nodes.PH(:,k1) = PHmat(:,strmatch(temp2{2},phlist,'exact'));
        end
    end
    
end
% Number of nodes
nnode = size(nodelinearray,2);
nodes.nnode = nnode;

% network.nnode = nodes.nnode;
% network.nphases = nodes.phases;
% network.NPH = nodes.PH;

%% Lines

% Iterate through lines in configuration file. Lines are added in order of
% their placement in the configuration file.
for k1 = 1:size(linelinearray,2)
    temp1 = strsplit(linelinearray{k1},' ');
    temp1 = temp1(2:end);
    
    for k2 = 1:length(temp1)    
        temp2 = strsplit(temp1{k2},'=');
        
        % Line TX (sending) node
        if strcmp('TXnode',temp2{1})
            % TX node name placed into string array
            TXnodename = temp2{2};
            lines.TXname{k1} = TXnodename;
            % TX node number (index) placed in array
            TXnodenum = strmatch(TXnodename,nodes.nodelist,'exact');
            lines.TXnum(k1) = TXnodenum;        
        end

        % Line RX (receiving) node
        if strcmp('RXnode',temp2{1})
            % RX node name placed into string array
            RXnodename = temp2{2};
            lines.RXname{k1} = RXnodename;
            % RX node number (index) placed into array
            RXnodenum = strmatch(RXnodename,nodes.nodelist,'exact');
            lines.RXnum(k1) = RXnodenum;        
        end

        % Line phase list and matrix
        if strcmp('phases',temp2{1})
            lines.phases{k1} = temp2{2};
            lines.PH(:,k1) = PHmat(:,strmatch(temp2{2},phlist,'exact'));
        end

        % Line configuration array
        if strcmp('config',temp2{1})
            lines.config{k1} = temp2{2}; 
        end

        % Line length [m]
        if strcmp('length',temp2{1})
            lines.length(k1) = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('unit',temp2{1})
                if strcmp('km',temp2{2})
                    lines.length(k1) = lines.length(k1)*1000;
                end
                if strcmp('ft',temp2{2})
                    lines.length(k1) = lines.length(k1)*0.3048;
                end
                if strcmp('mi',temp2{2})
                    lines.length(k1) = lines.length(k1)*1609.34;
                end
            end
        end
        
    end
    
    % Node connection matrix
    nodes.FM(TXnodenum,RXnodenum) = 1;
    nodes.FM(RXnodenum,TXnodenum) = -1;
        
end
% Number of lines
nline = size(linelinearray,2);
lines.nline = nline;

% network.nline = lines.nline;
% network.TXname = lines.TXname;
% network.TXnum = lines.TXnum;
% network.RXname = lines.RXname;
% network.RXnum = lines.RXnum;
% network.lphases = lines.phases;
% network.LPH = lines.PH;
% network.linelength = lines.length;

%% Find lines entering and leaving nodes

% Compute number of lines entering and leaving each node
intemp = zeros(1,nnode);
outtemp = zeros(1,nnode);
for k1 = 1:nline
    intemp(1,lines.RXnum(k1)) = intemp(1,lines.RXnum(k1)) + 1;
    outtemp(1,lines.TXnum(k1)) = outtemp(1,lines.TXnum(k1)) + 1;
end

% Matrix of lines coming into nodes, with nodes as columns, and lines as
% rows. Entries that are 0 are nonexistant (padded)
nodes.inlines = zeros(1,nnode);
nodes.innodes = zeros(1,nnode);
% Matrix of lines leaving nodes, with nodes as columns, and lines as
% rows. Entries that are 0 are nonexistant (padded)
nodes.outlines = zeros(1,nnode);
nodes.outnodes = zeros(1,nnode);
for k1 = 1:nnode
    incount = 1;
    outcount = 1;
    for k2 = 1:nline
        if k1 == lines.RXnum(k2)
            nodes.inlines(incount,k1) = k2;
            nodes.innodes(incount,k1) = lines.TXnum(k2);
            incount = incount + 1;
        end
        if k1 == lines.TXnum(k2)
            nodes.outlines(outcount,k1) = k2;
            nodes.outnodes(outcount,k1) = lines.RXnum(k2);
            outcount = outcount + 1;
        end
    end    
end

%% Line configurations

% Impedance matrices for configs [pu/m]
configs.FZpupl = zeros(3,3,size(configlinearray,2));

% Iterate through line configurations
for k1 = 1:size(configlinearray,2)
    temp1 = strsplit(configlinearray{k1},' ');
    temp1 = temp1(2:end);
    
    % Impedance (resistance and reactance) values for current line configuration
    for k2 = 1:length(temp1)    
        temp2 = strsplit(temp1{k2},'=');
        
        if strcmp('config',temp2{1})
            configs.conflist{k1} = temp2{2};
        end

        if strcmp('raa',temp2{1})
            raa = str2double(temp2{2});
        end

        if strcmp('xaa',temp2{1})
            xaa = str2double(temp2{2});
        end

        if strcmp('rab',temp2{1})
            rab = str2double(temp2{2});
        end

        if strcmp('xab',temp2{1})
            xab = str2double(temp2{2});
        end

        if strcmp('rac',temp2{1})
            rac = str2double(temp2{2});
        end

        if strcmp('xac',temp2{1})
            xac = str2double(temp2{2});
        end

        if strcmp('rbb',temp2{1})
            rbb = str2double(temp2{2});
        end

        if strcmp('xbb',temp2{1})
            xbb = str2double(temp2{2});
        end

        if strcmp('rbc',temp2{1})
            rbc = str2double(temp2{2});
        end

        if strcmp('xbc',temp2{1})
            xbc = str2double(temp2{2});
        end

        if strcmp('rcc',temp2{1})
            rcc = str2double(temp2{2});
        end

        if strcmp('xcc',temp2{1})
            xcc = str2double(temp2{2});
        end
        
        % Impedance per unit length multiplying factor for given unit
        if strcmp('unit',temp2{1}) || strcmp('lengthunit',temp2{1})
            if strcmp('pu/m',temp2{2})
                FZmult = 1;
            end
            if strcmp('pu/km',temp2{2})
                FZmult = 1e-3;
            end
            if strcmp('pu/ft',temp2{2})
                FZmult = 0.3048;
            end
            if strcmp('pu/mile',temp2{2})
                FZmult = 1/1609.34;
            end
            if strcmp('ohm/m',temp2{2})
                FZmult = 1/Zbase;
            end
            if strcmp('ohm/km',temp2{2})
                FZmult = 1e-3/Zbase;
            end
            if strcmp('ohm/ft',temp2{2})
                FZmult = 0.3048/Zbase;
            end
            if strcmp('ohm/mile',temp2{2})
                FZmult = 1/1609.34/Zbase;
            end
        end
                
    end
    
    % 3x3 impedance per unit length matrix for current line config [pu/m]
    configs.FZpupl(:,:,k1) = FZmult*...
        [(raa + 1j*xaa) (rab + 1j*xab) (rac + 1j*xac);
        (rab + 1j*xab) (rbb + 1j*xbb) (rbc + 1j*xbc);
        (rac + 1j*xac) (rbc + 1j*xbc) (rcc + 1j*xcc)];
    
end

% network.conflist = configs.conflist;
% network.FZpl = configs.FZpl;

clear FZmult

%% Reconcile lines and configs

% Impednace matrices for line [pu]
lines.FZpu = zeros(3,3,nline);

% Iterate through lines
for k1 = 1:nline
    
    % Find configuration for current line
    confnum = find(strcmp(lines.config{k1},configs.conflist));
    % Multiply 3x3 per unit length impedance matrix [ohm/m] by line length [m]
    % to obtain 3x3 line impdance matrix [pu]
    lines.FZpu(:,:,k1) = lines.length(k1)*configs.FZpupl(:,:,confnum);
    % Zero out impedance matrix values for phases that do not exist on line
    if lines.PH(1,k1) == 0
        lines.FZpu(1,:,k1) = 0; lines.FZpu(:,1,k1) = 0;
    end
    if lines.PH(2,k1) == 0
        lines.FZpu(2,:,k1) = 0; lines.FZpu(:,2,k1) = 0;
    end
    if lines.PH(3,k1) == 0
        lines.FZpu(3,:,k1) = 0; lines.FZpu(:,3,k1) = 0;
    end
          
end

% Resistance, and reactance matrices for all lines [pu]
lines.FRpu = real(lines.FZpu);
lines.FXpu = imag(lines.FZpu);

% network.FZpu = lines.FZpu;
% network.FRpu = lines.FRpu;
% network.FXpu = lines.FXpu;

%% Calculate Admittance matrices

% Compute per unit admittance matrices for all lines
% Admittance matrices for lines [pu]
lines.FYpu = zeros(3,3,nline);
for k1 = 1:nline
    if strcmp(lines.phases{k1},'a')
        lines.FYpu(1,1,k1) = 1./lines.FZpu(1,1,k1);
    elseif strcmp(lines.phases{k1},'b')
        lines.FYpu(2,2,k1) = 1./lines.FZpu(2,2,k1);
    elseif strcmp(lines.phases{k1},'c')
        lines.FYpu(3,3,k1) = 1./lines.FZpu(3,3,k1);
    elseif strcmp(lines.phases{k1},'ab')
        lines.FYpu(1:2,1:2,k1) = pinv(lines.FZpu(1:2,1:2,k1));
    elseif strcmp(lines.phases{k1},'bc')
    	lines.FYpu(2:3,2:3,k1) = pinv(lines.FZpu(2:3,2:3,k1));
    elseif strcmp(lines.phases{k1},'ac')
        temp = lines.FZpu(:,:,k1);
        temp(2,:) = []; temp(:,2) = [];
        temp = pinv(temp);
        lines.FYpu(1,1,k1) = temp(1,1);
        lines.FYpu(1,3,k1) = temp(1,2);
        lines.FYpu(3,1,k1) = temp(2,1);
        lines.FYpu(3,3,k1) = temp(2,2); 
    elseif strcmp(lines.phases{k1},'abc')
        lines.FYpu(:,:,k1) = pinv(lines.FZpu(:,:,k1));
    end
end
% Conductance and susceptance matrices for lines [pu]
lines.FGpu = real(lines.FYpu);
lines.FBpu = imag(lines.FYpu);

% network.FYpu = lines.FYpu;
% network.FGpu = lines.FGpu;
% network.FBpu = lines.FBpu;

%% Loads

% Load parameters
% Load real component [pu]
loads.ppu = zeros(3,nnode);
% Load reactive component [pu]
loads.qpu = zeros(3,nnode);
% Load constant power coefficient
loads.aPQ = zeros(3,nnode);
% Load constant current coefficient
loads.aI = zeros(3,nnode);
% Load constant impedance coefficient
loads.aZ = zeros(3,nnode);

% Iterate through loads in configuration
for k1 = 1:size(loadlinearray,2)
    temp1 = strsplit(loadlinearray{k1},' ');
    temp1 = temp1(2:end);
    
    for k2 = 1:length(temp1)    
        temp2 = strsplit(temp1{k2},'=');
    
        % Load node
        if strcmp('nodename',temp2{1})
            knode = strmatch(temp2{2},nodes.nodelist,'exact');
        end        
        if strcmp('nodenum',temp2{1})
            knode = str2double(temp2{2});
        end

        % Load connection - not used
        if strcmp('conn',temp2{1})
            if strmatch('wye',temp2{2})

            end
            if strmatch('delta',temp2{2})

            end
        end

        % Load phase(s)
        if strcmp('phases',temp2{1})
            kph = find(strcmp(temp2{2},{'a','b','c'}));
        end

        % Load type - not used
        if strcmp('type',temp2{1})
            loads.type{kph,knode} = temp2{2};
        end

        % Load constant power coefficient
        if strcmp('apq',temp2{1})
            loads.aPQ(kph,knode) = str2double(temp2{2});
        end

        % Load constant current coefficient
        if strcmp('ai',temp2{1})
            loads.aI(kph,knode) = str2double(temp2{2});
        end

        % Load constant impedance coefficient
        if strcmp('az',temp2{1})
            loads.aZ(kph,knode) = str2double(temp2{2});
        end

        % Load real component [pu]
        if strcmp('real',temp2{1})
            loads.ppu(kph,knode) = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('pu',temp2{2})
                loads.ppu(kph,knode) = loads.ppu(kph,knode);
            end
            if strcmp('W',temp2{2})
                loads.ppu(kph,knode) = loads.ppu(kph,knode)*1/Sbase;
            end
            if strcmp('kW',temp2{2})
                loads.ppu(kph,knode) = loads.ppu(kph,knode)*1e3/Sbase;
            end
            if strcmp('MW',temp2{2})
                loads.ppu(kph,knode) = loads.ppu(kph,knode)*1e6/Sbase;
            end
        end

        % Load reactive component [pu]
        if strcmp('reac',temp2{1})
            loads.qpu(kph,knode) = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('pu',temp2{2})
                loads.qpu(kph,knode) = loads.qpu(kph,knode);
            end
            if strcmp('VAr',temp2{2})
                loads.qpu(kph,knode) = loads.qpu(kph,knode)*1e0/Sbase;
            end
            if strcmp('kVAr',temp2{2})
                loads.qpu(kph,knode) = loads.qpu(kph,knode)*1e3/Sbase;
            end
            if strcmp('MVAr',temp2{2})
                loads.qpu(kph,knode) = loads.qpu(kph,knode)*1e6/Sbase;
            end
        end
        
    end
    
end
% Complex loads [pu]
loads.spu = loads.ppu + 1j*loads.qpu;

% network.loadtype = loads.type;
% network.loadconn = loads.conn;
% network.spu = loads.spu;
% network.ppu = loads.ppu;
% network.qpu = loads.qpu;
% network.aPQ = loads.aPQ;
% network.aI = loads.aI;
% network.aZ = loads.aZ;

%% Capacitors

% Capacitor parameters
% Capacitor matrix [pu]
caps.cappu = zeros(3,nnode);

% Iterate through capacitors in configuration file
for k1 = 1:size(capacitorlinearray,2)
    temp1 = strsplit(capacitorlinearray{k1},' ');
    temp1 = temp1(2:end);
    
    for k2 = 1:length(temp1)    
        temp2 = strsplit(temp1{k2},'=');
    
        % Capacitor node
        if strcmp('nodename',temp2{1})
            knode = strmatch(temp2{2},nodes.nodelist);
        end
        if strcmp('nodenum',temp2{1})
            knode = str2double(temp2{2});
        end

        % Capacitor connection - not used
        if strcmp('conn',temp2{1})
            if strcmp('wye',temp2{2})

            end
        end

        % Capacitor phase(s)
        if strcmp('phase',temp2{1})
            kph = strmatch(temp2{2},{'a','b','c'},'exact');
        end

        % Capacitance [pu]
        if strcmp('reac',temp2{1})
            caps.cappu(kph,knode) = str2double(temp2{2});
            temp2 = strsplit(temp1{5},'=');
            if strcmp('pu',temp2{2})
                caps.cappu(kph,knode) = caps.cappu(kph,knode);
            end
            if strcmp('VAr',temp2{2})
                caps.cappu(kph,knode) = caps.cappu(kph,knode)*1e0/Sbase;
            end
            if strcmp('kVAr',temp2{2})
                caps.cappu(kph,knode) = caps.cappu(kph,knode)*1e3/Sbase;
            end
            if strcmp('MVAr',temp2{2})
                caps.cappu(kph,knode) = caps.cappu(kph,knode)*1e6/Sbase;
            end
        end
    end   
    
end

% network.capstype = caps.type
% network.capsconn = caps.conn
% network.cappu = caps.cappu;

%% Controllers (cons)

% Controller parameters
% Controller apparent power capacity [pu]
cons.wmaxpu = zeros(3,nnode);
% ES controller frequency
cons.fes = zeros(3,nnode);
% ES controller high pass filter (hpf) frequency
cons.hpfes = zeros(3,nnode);
% ES controller low pass filter (lpf) frequency
cons.lpfes = zeros(3,nnode);
% ES controller integrator gain
cons.kintes = zeros(3,nnode);

% Iterate through cons in configuration file
for k1 = 1:size(controllerlinearray,2)
    temp1 = strsplit(controllerlinearray{k1},' ');
    temp1 = temp1(2:end);
    
    for k2 = 1:length(temp1)    
        temp2 = strsplit(temp1{k2},'=');

        % Controller node
        if strcmp('nodename',temp2{1})
            knode = find(strcmp(temp2{2},nodes.nodelist));
        end
        if strcmp('nodenum',temp2{1})
            knode = str2double(temp2{2});
        end

        % Controller connection - not used
        if strcmp('conn',temp2{1})
            if strcmp('wye',temp2{2})

            end
        end

        % Controller phases(s)
        if strcmp('phases',temp2{1})
            kph = find(strcmp(temp2{2},{'a','b','c'}));
        end

        % Controller apparent power capacity [pu]
        if strcmp('mag',temp2{1})
            cons.wmaxpu(kph,knode) = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('pu',temp2{2})
                cons.wmaxpu(kph,knode) = cons.wmaxpu(kph,knode);
            end
            if strcmp('VAr',temp2{2})
                cons.wmaxpu(kph,knode) = cons.wmaxpu(kph,knode)*1e0/Sbase;
            end
            if strcmp('kVAr',temp2{2})
                cons.wmaxpu(kph,knode) = cons.wmaxpu(kph,knode)*1e3/Sbase;
            end
            if strcmp('MVAr',temp2{2})
                cons.wmaxpu(kph,knode) = cons.wmaxpu(kph,knode)*1e6/Sbase;
            end
        end

        % ES controller frequency
        if strcmp('fes',temp2{1})
            cons.fes(kph,knode) = str2double(temp2{2});
        end

        % ES controller high pass filter (hpf) frequency
        if strcmp('hpfes',temp2{1})
            cons.hpfes(kph,knode) = str2double(temp2{2});
        end

        % ES controller low pass filter (lpf) frequency
        if strcmp('lpfes',temp2{1})
            cons.lpfes(kph,knode) = str2double(temp2{2});
        end

        % ES controller integrator gain
        if strcmp('kintes',temp2{1})
            cons.kintes(kph,knode) = str2double(temp2{2});
        end
        
    end
    
end

% Controller dispatch [pu]
cons.wpu = zeros(3,nnode);

% network.constype = cons.type;
% network.consconn = cons.conn;
% network.wmaxpu = cons.wmaxpu;
% network.fes = cons.fes;
% network.hpfes = cons.hpfes;
% network.lpfes = cons.lpfes;
% network.kintes = cons.kintes;

%% VVC

vvc.state = zeros(3,nnode);
vvc.type = zeros(3,nnode);
% VVC minimum voltage [pu]
vvc.Vminpu = zeros(3,nnode);
% VVC maximum voltage [pu]
vvc.Vmaxpu = zeros(3,nnode);
% VVC maximum reactive power [pu]
vvc.qminpu = zeros(3,nnode);
% VVC maximum reactive power [pu]
vvc.qmaxpu = zeros(3,nnode);

% Iterate through vvc in configuration file
for k1 = 1:size(vvclinearray,2)
    temp1 = strsplit(vvclinearray{k1},' ');
    temp1 = temp1(2:end);
    
    for k2 = 1:length(temp1)    
        temp2 = strsplit(temp1{k2},'=');
        
        % VVC node
        if strcmp('nodename',temp2{1})
            knode = find(strcmp(temp2{2},nodes.nodelist));
        end
        if strcmp('nodenum',temp2{1})
            knode = str2double(temp2{2});
        end

        % VVC connection - not used
        if strcmp('conn',temp2{1})
            if strcmp('wye',temp2{2})

            end
        end

        % VVC phases(s)
        if strcmp('phases',temp2{1})
            kph = find(strcmp(temp2{2},{'a','b','c'}));
        end
        
        % VVC state
        if strcmp('state',temp2{1})
            vvc.state(kph,knode) = str2double(temp2{2});
        end
        
        % VVC type
        if strcmp('type',temp2{1})
            vvc.type(kph,knode) = str2double(temp2{2});
        end
        
        % VVC minimum voltage [pu]
        if strcmp('Vmin',temp2{1})
            vvc.Vminpu(kph,knode) = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('pu',temp2{2})
                vvc.Vminpu(kph,knode) = vvc.Vminpu(kph,knode);
            end
            if strcmp('V',temp2{2})
                vvc.Vminpu(kph,knode) = vvc.Vminpu(kph,knode)*1e0/Vbase;
            end
            if strcmp('kV',temp2{2})
                vvc.Vminpu(kph,knode) = vvc.Vminpu(kph,knode)*1e3/Vbase;
            end
            if strcmp('MV',temp2{2})
                vvc.Vminpu(kph,knode) = vvc.Vminpu(kph,knode)*1e6/Vbase;
            end
        end
        
        % VVC maximum voltage [pu]
        if strcmp('Vmax',temp2{1})
            vvc.Vmaxpu(kph,knode) = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('pu',temp2{2})
                vvc.Vmaxpu(kph,knode) = vvc.Vmaxpu(kph,knode);
            end
            if strcmp('V',temp2{2})
                vvc.Vmaxpu(kph,knode) = vvc.Vmaxpu(kph,knode)*1e0/Vbase;
            end
            if strcmp('kV',temp2{2})
                vvc.Vmaxpu(kph,knode) = vvc.Vmaxpu(kph,knode)*1e3/Vbase;
            end
            if strcmp('MV',temp2{2})
                vvc.Vmaxpu(kph,knode) = vvc.Vmaxpu(kph,knode)*1e6/Vbase;
            end
        end
                
        % VVC minimum reactive power [pu]
        if strcmp('qmin',temp2{1})
            vvc.qminpu(kph,knode) = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('pu',temp2{2})
                vvc.qminpu(kph,knode) = vvc.qminpu(kph,knode);
            end
            if strcmp('VAr',temp2{2})
                vvc.qminpu(kph,knode) = vvc.qminpu(kph,knode)*1e0/Sbase;
            end
            if strcmp('kVAr',temp2{2})
                vvc.qminpu(kph,knode) = vvc.qminpu(kph,knode)*1e3/Sbase;
            end
            if strcmp('MVAr',temp2{2})
                vvc.qminpu(kph,knode) = vvc.qminpu(kph,knode)*1e6/Sbase;
            end
        end
        
        % VVC maximum reactive power [pu]
        if strcmp('qmin',temp2{1})
            vvc.qmaxpu(kph,knode) = str2double(temp2{2});
            temp2 = strsplit(temp1{k2+1},'=');
            if strcmp('pu',temp2{2})
                vvc.qmaxpu(kph,knode) = vvc.qmaxpu(kph,knode);
            end
            if strcmp('VAr',temp2{2})
                vvc.qmaxpu(kph,knode) = vvc.qmaxpu(kph,knode)*1e0/Sbase;
            end
            if strcmp('kVAr',temp2{2})
                vvc.qmaxpu(kph,knode) = vvc.qmaxpu(kph,knode)*1e3/Sbase;
            end
            if strcmp('MVAr',temp2{2})
                vvc.qmaxpu(kph,knode) = vvc.qmaxpu(kph,knode)*1e6/Sbase;
            end
        end
        
    end
end

% VVC dispatch [pu]
vvc.vvcpu = zeros(3,nnode);

%% Place all structs into network struct

network.base = base;
network.nodes = nodes;
network.lines = lines;
network.configs = configs;
network.loads = loads;
network.caps = caps;
network.cons = cons;
network.vvc = vvc;

%%

clear linearray baselinearray nodelinearray linelinearray configlinearray
clear loadlinearray capacitorlinearray controllerlinearray vvclinearray
clear templinearray
clear kph knode k1 ph num tline ans PHmat
clear phlist confnum fid fn fp intemp k1 k2 temp1 temp2 outtemp
clear raa xaa rab xab rac xac rbb xbb rbc xbc rcc xcc
clear TXnodename TXnodenum RXnodename RXnodenum intemp outtemp
clear incount inlines innodes outcount outlines outnodes
