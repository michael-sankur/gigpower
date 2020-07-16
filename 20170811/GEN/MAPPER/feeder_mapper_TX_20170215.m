% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

clc, clear all, close all

%% Load all feeder data

feeder.name = '2node_fullphase';

% Read from csv file, input first 4 columns as strings (general case)
% Change file path as needed
fp = ['C:\Users\Michael\Desktop\reconfiguration\Feeders\'];
fn = 'ieee_13node_new.txt';
fn = 'ieee_37node_singlephase.txt';

linearray = {};

fid = fopen([fp fn]);
tline = fgetl(fid);
while ischar(tline)
    linearray{end+1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

baselinearray = {};
nodelinearray = {};
linelinearray = {};
configlinearray = {};
loadlinearray = {};
capacitorlinearray = {};
controllerlinearray = {};
for k1 = 1:size(linearray,2)
    templinearray = strsplit(linearray{k1},' ');
    if strmatch('base', templinearray{1},'exact')
        baselinearray{end+1} = linearray{k1};
    end
    if strmatch('node', templinearray{1},'exact')
        nodelinearray{end+1} = linearray{k1};
    end
    if strmatch('line', templinearray{1},'exact')
        linelinearray{end+1} = linearray{k1};
    end
    if strmatch('config', templinearray{1},'exact')
        configlinearray{end+1} = linearray{k1};
    end
    if strmatch('load', templinearray{1},'exact')
        loadlinearray{end+1} = linearray{k1};
    end
    if strmatch('capacitor', templinearray{1},'exact')
        capacitorlinearray{end+1} = linearray{k1};
    end
    if strmatch('controller', templinearray{1},'exact')
        controllerlinearray{end+1} = linearray{k1};
    end
end

%% Base Values

for k1 = 1:size(baselinearray,2)
    temp1 = strsplit(baselinearray{k1},' ');
    temp1 = temp1(2:end);
    
    temp2 = strsplit(temp1{1},'=');
    if strmatch('Vbase',temp2{1},'exact')
        feeder.Vbase = str2double(temp2{2});
        temp2 = strsplit(temp1{2},'=');
        if strmatch('V',temp2{2},'exact')
            feeder.Vbase = feeder.Vbase;
        end
        if strmatch('kV',temp2{2},'exact')
            feeder.Vbase = 1000*feeder.Vbase;
        end
        if strmatch('MV',temp2{2},'exact')
            feeder.Vbase = 1e6*feeder.Vbase;
        end
    end
    
    if strmatch('Sbase',temp2{1},'exact')
        feeder.Sbase = str2double(temp2{2});
        temp2 = strsplit(temp1{2},'=');
        if strmatch('VAr',temp2{2},'exact')
            feeder.Sbase = feeder.Sbase;
        end
        if strmatch('kVAr',temp2{2},'exact')
            feeder.Sbase = 1000*feeder.Sbase;
        end
        if strmatch('MVAr',temp2{2},'exact')
            feeder.Sbase = 1e6*feeder.Sbase;
        end
    end
    
end

feeder.Ibase = feeder.Sbase/feeder.Vbase; % [A]
feeder.Zbase = feeder.Vbase/feeder.Ibase; % [Ohm]

%% Nodes

phlist = {'a','b','c','ab','bc','ac','abc'};
PHmat = [1 0 0;
    0 1 0;
    0 0 1;
    1 1 0;
    0 1 1;
    1 0 1;
    1 1 1]';

for k1 = 1:size(nodelinearray,2)
    temp1 = strsplit(nodelinearray{k1},' ');
    temp1 = temp1(2:end);
        
    temp2 = strsplit(temp1{1},'=');
    if strmatch('nodename',temp2{1},'exact')
        nodes.nodelist(k1) = temp2(2);
    end 
    
    temp2 = strsplit(temp1{2},'=');
    if strmatch('phases',temp2{1},'exact')
        nodes.phases{k1} = temp2{2};
        nodes.PH(:,k1) = PHmat(:,strmatch(temp2{2},phlist,'exact'));
    end
    
end
nodes.nnode = k1;

%% Lines

% lines.tx = []; lines.TX = [];
% lines.rx = []; lines.RX = [];
for k1 = 1:size(linelinearray,2)
    temp1 = strsplit(linelinearray{k1},' ');
    temp1 = temp1(2:end);
    
    temp2 = strsplit(temp1{1},'=');
    if strmatch('TXnode',temp2{1},'exact')
        TXnodename = temp2{2};
        TXnodenum = strmatch(TXnodename,nodes.nodelist,'exact');
        lines.TXnum(k1) = TXnodenum;
        lines.TXname{k1} = TXnodename;
    end
    
    temp2 = strsplit(temp1{2},'=');
    if strmatch('RXnode',temp2{1},'exact')
        RXnodename = temp2{2};
        RXnodenum = strmatch(RXnodename,nodes.nodelist,'exact');
        lines.RXnum(k1) = RXnodenum;
        lines.RXname{k1} = RXnodename;
    end
    
    temp2 = strsplit(temp1{3},'=');
    if strmatch('phases',temp2{1},'exact')
        lines.phases{k1} = temp2{2};
        lines.PH(:,k1) = PHmat(:,strmatch(temp2{2},phlist,'exact'));
    end
    
    temp2 = strsplit(temp1{4},'=');
    if strmatch('config',temp2{1},'exact')
        lines.config{k1} = temp2{2}; 
    end
    
    temp2 = strsplit(temp1{5},'=');
    if strmatch('length',temp2{1},'exact')
        lines.length(k1) = str2double(temp2{2}); 
    end
    
    temp2 = strsplit(temp1{6},'=');
    if strmatch('unit',temp2{1},'exact')
        if strmatch('km',temp2{2},'exact')
            lines.length(k1) = lines.length(k1)/1000;
        end
        if strmatch('ft',temp2{2},'exact')
            lines.length(k1) = lines.length(k1)*0.3048;
        end
        if strmatch('mi',temp2{2},'exact')
            lines.length(k1) = lines.length(k1)/1609.34;
        end
    end
    
    nodes.FM(TXnodenum,RXnodenum) = 1;
    nodes.FM(RXnodenum,TXnodenum) = -1;
        
end
lines.nline = k1;

inmat = [];
outmat = [];

for k1 = 1:nodes.nnode    
    
    intemp = [];
    outtemp = [];
   
    for k2 = 1:lines.nline
        
        if k1 == lines.RXnum(k2)
            intemp(end+1,1) = k2;
        end
        if k1 == lines.TXnum(k2)
            outtemp(end+1,1) = k2;
        end        
        
    end
    
    if isempty(intemp)
        intemp = 0;
    end    
    if size(intemp,1) < size(inmat,1)
        intemp = [intemp; zeros(size(inmat,1)-size(intemp,1),1)];
        inmat(:,k1) = intemp;
    elseif size(intemp,1) == size(inmat,1)
        inmat(:,k1) = intemp;
    elseif size(intemp,1) > size(inmat,1)
        inmat = [inmat; zeros(size(intemp,1)-size(inmat,1),size(inmat,2))];
        inmat(:,k1) = intemp;
    end
    
    if isempty(outtemp)
        outtemp = 0;
    end    
    if size(outtemp,1) < size(outmat,1)
        outtemp = [outtemp; zeros(size(outmat,1)-size(outtemp,1),1)];
        outmat(:,k1) = outtemp;
    elseif size(outtemp,1) == size(outmat,1)
        outmat(:,k1) = outtemp;
    elseif size(outtemp,1) > size(outmat,1)
        outmat = [outmat; zeros(size(outtemp,1)-size(outmat,1),size(outmat,2))];
        outmat(:,k1) = outtemp;
    end
    
end

nodes.inmat = inmat;
nodes.outmat = outmat;

clear TXnodename TXnodenum RXnodename RXnodenum inmat outmat

%% Line config

for k1 = 1:size(configlinearray,2)
    temp1 = strsplit(configlinearray{k1},' ');
    temp1 = temp1(2:end);
    
    temp2 = strsplit(temp1{1},'=');
    if strmatch('config',temp2{1},'exact')
        configs.conflist{k1} = temp2{2};
    end
    
    temp2 = strsplit(temp1{3},'=');
    if strmatch('raa',temp2{1},'exact')
        raa = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{4},'=');
    if strmatch('xaa',temp2{1},'exact')
        xaa = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{5},'=');
    if strmatch('rab',temp2{1},'exact')
        rab = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{6},'=');
    if strmatch('xab',temp2{1},'exact')
        xab = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{7},'=');
    if strmatch('rac',temp2{1},'exact')
        rac = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{8},'=');
    if strmatch('xac',temp2{1},'exact')
        xac = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{9},'=');
    if strmatch('rbb',temp2{1},'exact')
        rbb = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{10},'=');
    if strmatch('xbb',temp2{1},'exact')
        xbb = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{11},'=');
    if strmatch('rbc',temp2{1},'exact')
        rbc = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{12},'=');
    if strmatch('xbc',temp2{1},'exact')
        xbc = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{13},'=');
    if strmatch('rcc',temp2{1},'exact')
        rcc = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{14},'=');
    if strmatch('xcc',temp2{1},'exact')
        xcc = str2double(temp2{2});
    end
    
    configs.FZpl(:,:,k1) = [(raa + 1j*xaa) (rab + 1j*xab) (rac + 1j*xac);
        (rab + 1j*xab) (rbb + 1j*xbb) (rbc + 1j*xbc);
        (rac + 1j*xac) (rbc + 1j*xbc) (rcc + 1j*xcc)];
            
    temp2 = strsplit(temp1{2},'=');
    if strmatch('unit',temp2{1},'exact')
        if strmatch('ohm/m',temp2{2},'exact')
            configs.FZpl(:,:,k1) = configs.FZpl(:,:,k1);
        end
        if strmatch('ohm/km',temp2{2},'exact')
            configs.FZpl(:,:,k1) = configs.FZpl(:,:,k1)/1000;
        end
        if strmatch('ohm/ft',temp2{2},'exact')
            configs.FZpl(:,:,k1) = configs.FZpl(:,:,k1)*0.3048;
        end
        if strmatch('ohm/mile',temp2{2},'exact')
            configs.FZpl(:,:,k1) = configs.FZpl(:,:,k1)/1609.34;
        end
    end
    
end

clear raa xaa rab xab rac xac rbb xbb rbc xbc rcc xcc

%% Reconcile lines and configs

phlist = {'a','b','c','ab','bc','ac','abc'};

for k1 = 1:lines.nline
    
    confnum = strmatch(lines.config{k1},configs.conflist,'exact');
    lines.FZ(:,:,k1) = lines.length(k1)*configs.FZpl(:,:,confnum);
    if lines.PH(1,k1) == 0
        lines.FZ(1,:,k1) = 0; lines.FZ(:,1,k1) = 0;
    end
    if lines.PH(2,k1) == 0
        lines.FZ(2,:,k1) = 0; lines.FZ(:,2,k1) = 0;
    end
    if lines.PH(3,k1) == 0
        lines.FZ(3,:,k1) = 0; lines.FZ(:,3,k1) = 0;
    end
    
    lines.FZpu(:,:,k1) = lines.FZ(:,:,k1)/feeder.Zbase;
        
end

lines.FYpu = zeros(3,3,lines.nline);
for k1 = 1:lines.nline
    if lines.PH(1,k1) == 1 && lines.PH(2,k1) == 0 && lines.PH(3,k1) == 0
        lines.FYpu(1,1,k1) = 1/lines.FZpu(1,1,k1);
    elseif lines.PH(1,k1) == 0 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 0
        lines.FYpu(2,2,k1) = 1./lines.FZpu(2,2,k1);
    elseif lines.PH(1,k1) == 0 && lines.PH(2,k1) == 0 && lines.PH(3,k1) == 1
        lines.FYpu(3,3,k1) = 1./lines.FZpu(3,3,k1);
    elseif lines.PH(1,k1) == 1 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 0
        lines.FYpu(1:2,1:2,k1) = inv(lines.FZpu(1:2,1:2,k1));
    elseif lines.PH(1,k1) == 0 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 1
    	lines.FYpu(2:3,2:3,k1) = inv(lines.FZpu(2:3,2:3,k1));
    elseif lines.PH(1,k1) == 1 && lines.PH(2,k1) == 0 && lines.PH(3,k1) == 1
        temp = lines.FZpu(:,:,k1);
        temp(2,:) = []; temp(:,2) = [];
        temp = inv(temp);
        lines.FYpu(1,1,k1) = temp(1,1);
        lines.FYpu(1,3,k1) = temp(1,2);
        lines.FYpu(3,1,k1) = temp(2,1);
        lines.FYpu(3,3,k1) = temp(2,2); 
    elseif lines.PH(1,k1) == 1 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 1
        lines.FYpu(:,:,k1) = inv(lines.FZpu(:,:,k1));
    end
end

%% Loads

loads.p = zeros(3,nodes.nnode);
loads.q = zeros(3,nodes.nnode);
loads.aPQ = zeros(3,nodes.nnode);
loads.aI = zeros(3,nodes.nnode);
loads.aV = zeros(3,nodes.nnode);
for k1 = 1:size(loadlinearray,2)
    temp1 = strsplit(loadlinearray{k1},' ');
    temp1 = temp1(2:end);
    
    temp2 = strsplit(temp1{1},'=');
    if strmatch('nodename',temp2{1},'exact')
        knode = strmatch(temp2{2},nodes.nodelist,'exact');
    end
    if strmatch('nodenum',temp2{1},'exact')
        knode = str2double(temp2{2});
    end
        
    temp2 = strsplit(temp1{2},'=');
    if strmatch('conn',temp2{1},'exact')
        if strmatch('wye',temp2{2},'exact')
%             loads.conn
        end
    end
    
    temp2 = strsplit(temp1{3},'=');
    if strmatch('phases',temp2{1},'exact')
        ph = temp2{2};
        kph = strmatch(temp2{2},{'a','b','c'},'exact');
    end
    
    temp2 = strsplit(temp1{4},'=');
    if strmatch('type',temp2{1},'exact')
        loads.type{kph,knode} = temp2{2};
    end
    
    temp2 = strsplit(temp1{5},'=');
    if strmatch('apq',temp2{1},'exact')
        loads.aPQ(kph,knode) = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{6},'=');
    if strmatch('ai',temp2{1},'exact')
        loads.aI(kph,knode) = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{7},'=');
    if strmatch('az',temp2{1},'exact')
        loads.aZ(kph,knode) = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{8},'=');
    if strmatch('real',temp2{1},'exact')
        loads.p(kph,knode) = str2double(temp2{2});
        temp2 = strsplit(temp1{9},'=');
        if strmatch('kW',temp2{2},'exact')
            loads.p(kph,knode) = loads.p(kph,knode)*1000;
        end
    end
    
    temp2 = strsplit(temp1{10},'=');
    if strmatch('reac',temp2{1},'exact')
        loads.q(kph,knode) = str2double(temp2{2});
        temp2 = strsplit(temp1{11},'=');
        if strmatch('kVAr',temp2{2},'exact')
            loads.q(kph,knode) = loads.q(kph,knode)*1000;
        end
    end
    
    loads.s = loads.p + 1j*loads.q;
    
end
loads.spu = loads.s/feeder.Sbase;

%% Capacitors

caps.cap = zeros(3,nodes.nnode);
for k1 = 1:size(capacitorlinearray,2)
    temp1 = strsplit(capacitorlinearray{k1},' ');
    temp1 = temp1(2:end);
    
    temp2 = strsplit(temp1{1},'=');
    if strmatch('nodename',temp2{1},'exact')
        knode = strmatch(temp2{2},nodes.nodelist,'exact');
    end
    if strmatch('nodenum',temp2{1},'exact')
        knode = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{2},'=');
    if strmatch('conn',temp2{1},'exact')
        if strmatch('wye',temp2{2},'exact')
%             loads.conn
        end
    end
    
    temp2 = strsplit(temp1{3},'=');
    if strmatch('phase',temp2{1},'exact')
        ph = temp2{2};
        kph = strmatch(temp2{2},{'a','b','c'},'exact');
    end
    
    temp2 = strsplit(temp1{4},'=');
    if strmatch('reac',temp2{1},'exact');
        caps.cap(kph,knode) = str2double(temp2{2})
        temp2 = strsplit(temp1{5},'=');
        if strmatch('kVAr',temp2{2},'exact')
            caps.cap(kph,knode) = caps.cap(kph,knode)*1000;
        end
    end
    
end
caps.cappu = caps.cap/feeder.Sbase;

%% Controllers

controllers.wmax = zeros(3,nodes.nnode);
for k1 = 1:size(controllerlinearray,2)
    temp1 = strsplit(controllerlinearray{k1},' ');
    temp1 = temp1(2:end);

    temp2 = strsplit(temp1{1},'=');
    if strmatch('nodename',temp2{1},'exact')
        knode = strmatch(temp2{2},nodes.nodelist,'exact');
    end
    if strmatch('nodenum',temp2{1},'exact')
        knode = str2double(temp2{2});
    end
    
    temp2 = strsplit(temp1{2},'=');
    if strmatch('conn',temp2{1},'exact')
        if strmatch('wye',temp2{2},'exact')
%             loads.conn
        end
    end
    
    temp2 = strsplit(temp1{3},'=');
    if strmatch('phases',temp2{1},'exact')
        ph = temp2{2};
        kph = strmatch(temp2{2},{'a','b','c'},'exact');
    end
    
    temp2 = strsplit(temp1{4},'=');
    if strmatch('mag',temp2{1},'exact');
        controllers.wmax(kph,knode) = str2double(temp2{2});
        temp2 = strsplit(temp1{5},'=');
        if strmatch('kVAr',temp2{2},'exact')
            controllers.wmax(kph,knode) = controllers.wmax(kph,knode)*1000;
        end
    end
    
end

controllers.wmaxpu = controllers.wmax/feeder.Sbase;

%%

clear linearray baselinearray nodelinearray linelinearray configlinearray
clear loadlinearray capacitorlinearray controllerlinearray templinearray
clear kph knode k1 ph num tline ans PHmat
clear phlist confnum fid fn fp intemp k1 k2 temp1 temp2 outtemp
