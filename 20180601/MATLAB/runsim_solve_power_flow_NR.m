% Michael Sankur - msankur@lbl.gov
% 2018.06.01

clc, clear, close all

% Add directory and subdirectories to path
path(path,genpath('C:\Users\kathl\Desktop\LinDist3Flow\20180601\MATLAB\'));

%% Load network

% home = pwd; %returns a string for the path of the base folder in which the .m file is housed.
% if ispc %This should identify the appropriate pathseparator for the operating system so that it will run on a windows, mac, or linux machine.
%     pathseparator = '\';
% else
%     pathseparator = '/';
% end
% genpathstr=strcat(home); %path string to the GEN folder
% feederpathstr=strcat(home,pathseparator, 'Networks', pathseparator); %path string to the Networks folder
% path(path,genpath(genpathstr));

% File path and file name of network configuration file
fp = 'C:\Users\kathl\Desktop\LinDist3Flow\20180601\NETWORKS\';

fp2 = 'C:\Users\kathl\Desktop\LinDist3Flow\20180601\testpvnum0';
fn = '05node_singlephase_radial.txt'; %5 nodes, each node has 3 phases
%fn = '05node_fullphase_radial.txt';
%fn = 'ieee_13node_balance.txt';
%fn = '03node_fullphase_mesh_example.txt'
%fn = '03node_singlephase_mesh_example.txt'
%fn = '03node_fullphase_radial_example.txt'
%fn = 'ieee_13node_mesh_open.txt';
%fn = '03node_fullphase_radial_COPY.txt';
%fn = 'ieee_34node - Copy.txt';

[network1] = network_mapper_function(fp, fn);

%% Network paramaters

nnode = network1.nodes.nnode; % Number of nodes
nline = network1.lines.nline; % Number of lines

%% Node parameters



%% Line parameters

% network1.lines.FZpu

%% Load parameters

% Change ZIP parameters for all loads/demands
% network1.loads.aPQ = 0.75*ones(3,nnode).*(network1.loads.spu ~= 0);
% network1.loads.aI = 0.10*ones(3,nnode).*(network1.loads.spu ~= 0);
% network1.loads.aZ = 0.15*ones(3,nnode).*(network1.loads.spu ~= 0);

network1.loads.aPQ = 1.00*ones(3,nnode).*(network1.loads.spu ~= 0);
network1.loads.aI = 0.00*ones(3,nnode).*(network1.loads.spu ~= 0);
network1.loads.aZ = 0.00*ones(3,nnode).*(network1.loads.spu ~= 0);

network1.loads.spu = 2*network1.loads.spu;


%% Phases That Exist
phase_count = 0;
a = 0;
b = 0;
c = 0;
for phaseCount=1:length(network1.nodes.phases)
    curr = network1.nodes.phases{phaseCount};
    if contains(curr, 'a') & a == 0
        phase_count = phase_count + 1;
        a = 1;
    end
    if contains(curr, 'b') & b == 0
        b = 1;
        phase_count = phase_count + 1;
    end
    if contains(curr, 'c') & c == 0
        phase_count = phase_count+1;
        c = 1;
    
    end
end


%% Adjusting spu

csv_file_count = 1;
all_demands= cell(phase_count,nnode);

for ph = 1:3
    for kn = 1:nnode
        if (a == 1 & ph == 1) | (b == 1 & ph == 2) | (c == 1 & ph == 3)
        csvfile = readmatrix('C:\Users\kathl\Desktop\LinDist3Flow\20180601\testpvnum0\node_'+string(csv_file_count)+'_pv_0_minute_normalized.csv');
        demand = csvfile(:, 6);
        all_demands{ph, kn} = demand;
        csv_file_count = csv_file_count + 1;
        end
    end
end
phstr = {'a','b','c'};
%for kt = 1:length(csvfile) %time
for kt = 1:50
    network1.loads.spu_nom = network1.loads.spu
    for ph = 1:3 %no of phases
        for kn = 1:nnode %no of nodes
            if (a == 1 & ph == 1) | (b == 1 & ph == 2) | (c == 1 & ph == 3)       
                network1.loads.spu_nom(ph,kn) = ...
                    network1.loads.spu(ph,kn)*all_demands{ph,kn}(kt)%(ph,kn,kt)
            end
        end
    end
    
    network1.cons.wpu = zeros(3,nnode);
    network1.vvc.vvcpu = zeros(3,nnode);
    slacknode = 1;
    Vslack = 1*[1;
        1*exp(1j*-120*pi/180);
        1*exp(1j*120*pi/180)];
    [NRRES0, VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, iter] = NR3(network1,slacknode,Vslack,[],[],1e-9);


    fp = 'C:\Users\kathl\Desktop\LinDist3Flow\20180601\MATLAB\';
    fn = 'NR3-matlab.txt';
    fid = fopen(fn,'a');

    fprintf(fid,'VNR\n\n');
    for ph = 1:3
        fprintf(fid,[phstr{ph} '\n']);
        for k1 = 1:nnode
            fprintf(fid,'%1.6f %1.6f\n',real(VNR0(ph,k1)),imag(VNR0(ph,k1)));
        end
    end

    fprintf(fid,'\n\nINR\n\n');
    for ph = 1:3
        fprintf(fid,[phstr{ph} '\n']);
        for k1 = 1:nline
            fprintf(fid,'%1.6f %1.6f\n',real(INR0(ph,k1)),imag(INR0(ph,k1)));
        end
    end

    fprintf(fid,'\n\nSTXNR\n\n');
    for ph = 1:3
        fprintf(fid,[phstr{ph} '\n']);
        for k1 = 1:nline
            fprintf(fid,'%1.6f %1.6f\n',real(STXNR0(ph,k1)),imag(STXNR0(ph,k1)));
        end
    end

    fprintf(fid,'\n\nSRXNR\n\n');
    for ph = 1:3
        fprintf(fid,[phstr{ph} '\n']);
        for k1 = 1:nline
            fprintf(fid,'%1.6f %1.6f\n',real(SRXNR0(ph,k1)),imag(SRXNR0(ph,k1)));
        end
    end

    fprintf(fid,'\n\niNR\n\n');
    for ph = 1:3
        fprintf(fid,[phstr{ph} '\n']);
        for k1 = 1:nnode
            fprintf(fid,'%1.6f %1.6f\n',real(iNR0(ph,k1)),imag(iNR0(ph,k1)));
        end
    end

    fprintf(fid,'\n\nsNR\n\n');
    for ph = 1:3
        fprintf(fid,[phstr{ph} '\n']);
        for k1 = 1:nnode
            fprintf(fid,'%1.6f %1.6f\n',real(sNR0(ph,k1)),imag(sNR0(ph,k1)));
        end
    end

    fclose(fid);



end

%% testing sth
fopen

%% Capacitor parameters

% network1.caps.cappu = 1*network1.caps.cappu;

%% Controller parameters

% network1.cons.wmaxpu = 1*network1.cons.wmaxpu;

% zero out control dispatch
network1.cons.wpu = zeros(3,nnode);

%% VVC parameters

% zero out vvc dispatch
network1.vvc.vvcpu = zeros(3,nnode);

%% Simulation parameters

% set slack node and slack node voltage
slacknode = 1;
Vslack = 1*[1;
    1*exp(1j*-120*pi/180);
    1*exp(1j*120*pi/180)];

%% Call NR Algorthm to solve power flow

% Solve power flow
[NRRES0, VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, iter] = NR3(network1,slacknode,Vslack,[],[],1e-9)

%%

phstr = {'a','b','c'};

% VNR0 = round(1e6*VNR0)/1e6

fp = 'C:\Users\kathl\Desktop\LinDist3Flow\20180601\MATLAB\';
fn = ['NR3-' erase(fn,'.txt') '-matlab' '.txt'];

fid = fopen([fp fn],'w');

fprintf(fid,'VNR\n\n');
for ph = 1:3
    fprintf(fid,[phstr{ph} '\n']);
    for k1 = 1:nnode
        fprintf(fid,'%1.6f %1.6f\n',real(VNR0(ph,k1)),imag(VNR0(ph,k1)));
    end
end

fprintf(fid,'\n\nINR\n\n');
for ph = 1:3
    fprintf(fid,[phstr{ph} '\n']);
    for k1 = 1:nline
        fprintf(fid,'%1.6f %1.6f\n',real(INR0(ph,k1)),imag(INR0(ph,k1)));
    end
end

fprintf(fid,'\n\nSTXNR\n\n');
for ph = 1:3
    fprintf(fid,[phstr{ph} '\n']);
    for k1 = 1:nline
        fprintf(fid,'%1.6f %1.6f\n',real(STXNR0(ph,k1)),imag(STXNR0(ph,k1)));
    end
end

fprintf(fid,'\n\nSRXNR\n\n');
for ph = 1:3
    fprintf(fid,[phstr{ph} '\n']);
    for k1 = 1:nline
        fprintf(fid,'%1.6f %1.6f\n',real(SRXNR0(ph,k1)),imag(SRXNR0(ph,k1)));
    end
end

fprintf(fid,'\n\niNR\n\n');
for ph = 1:3
    fprintf(fid,[phstr{ph} '\n']);
    for k1 = 1:nnode
        fprintf(fid,'%1.6f %1.6f\n',real(iNR0(ph,k1)),imag(iNR0(ph,k1)));
    end
end

fprintf(fid,'\n\nsNR\n\n');
for ph = 1:3
    fprintf(fid,[phstr{ph} '\n']);
    for k1 = 1:nnode
        fprintf(fid,'%1.6f %1.6f\n',real(sNR0(ph,k1)),imag(sNR0(ph,k1)));
    end
end

fclose(fid);