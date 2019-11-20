% Michael Sankur - msankur@lbl.gov
% 2018.06.01

clc, clear, close all

% Add directory and subdirectories to path
path(path,genpath('/home/michael/Dropbox/Unbalanced LinDistflow/20180601/MATLAB/'));

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
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180601/NETWORKS/';
fn = '03node_fullphase_radial_example.txt';
fn = '03node_singlephase_radial_example.txt';
% fn = 'ieee_13node_mesh_open.txt';

[network1] = network_mapper_function(fp, fn);

fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180601/NETWORKS/';
fn = '03node_fullphase_mesh_example.txt';
fn = '03node_singlephase_mesh_example.txt';
% fn = 'ieee_13node_mesh_closed.txt';

[network2] = network_mapper_function(fp, fn);

%% Network paramaters

nnode1 = network1.nodes.nnode; % Number of nodes
nline1 = network1.lines.nline; % Number of lines

nnode2 = network1.nodes.nnode; % Number of nodes
nline2 = network1.lines.nline; % Number of lines

%% Node parameters



%% Line parameters

% network1.lines.FZpu = 1/10*network1.lines.FZpu;

% network2.lines.FZpu = 1/10*network2.lines.FZpu;

%% Load parameters

% Change ZIP parameters for all loads/demands
% network1.loads.aPQ = 0.75*ones(3,nnode).*(network1.loads.spu ~= 0);
% network1.loads.aI = 0.10*ones(3,nnode).*(network1.loads.spu ~= 0);
% network1.loads.aZ = 0.15*ones(3,nnode).*(network1.loads.spu ~= 0);

network1.loads.aPQ = 1.00*ones(3,nnode1).*(network1.loads.spu ~= 0);
network1.loads.aI = 0.00*ones(3,nnode1).*(network1.loads.spu ~= 0);
network1.loads.aZ = 0.00*ones(3,nnode1).*(network1.loads.spu ~= 0);

network1.loads.spu(:,4) = 1.5*network1.loads.spu(:,3);

% network1.loads.spu = 2*network1.loads.spu;

network2.loads.aPQ = 1.00*ones(3,nnode2).*(network2.loads.spu ~= 0);
network2.loads.aI = 0.00*ones(3,nnode2).*(network2.loads.spu ~= 0);
network2.loads.aZ = 0.00*ones(3,nnode2).*(network2.loads.spu ~= 0);

network2.loads.spu(:,4) = 1.5*network2.loads.spu(:,3);

% network2.loads.spu(:,16:28) = 2*network2.loads.spu(:,16:28);

% network2.loads.spu = 2*network2.loads.spu;

%% Capacitor parameters

% network1.caps.cappu = 1*network1.caps.cappu;

% network2.caps.cappu = 1*network2.caps.cappu;

%% Controller parameters

% network1.cons.wmaxpu = 1*network1.cons.wmaxpu;

% zero out control dispatch
network1.cons.wpu = zeros(3,nnode1);

network2.cons.wpu = zeros(3,nnode2);

%% VVC parameters

% zero out vvc dispatch
network1.vvc.vvcpu = zeros(3,nnode1);

network2.vvc.vvcpu = zeros(3,nnode2);

%% Simulation parameters

% set slack node and slack node voltage
slacknode = 1;
Vslack = 1*[1;
    1*exp(1j*-120*pi/180);
    1*exp(1j*120*pi/180)];

%% Call NR Algorthm to solve power flow

% Solve power flow
[NRRES1, VNR1, INR1, STXNR1, SRXNR1, iNR1, sNR1, iter1] = NR3(network1,slacknode,Vslack,[],[],1e-9)

[NRRES2, VNR2, INR2, STXNR2, SRXNR2, iNR2, sNR2, iter2] = NR3(network2,slacknode,Vslack,[],[],1e-9)

%%

VNR1(:,3).*conj(network2.lines.FYpu(:,:,4)*(VNR1(:,3) - VNR1(:,4)))

%%

% phstr = {'a','b','c'};

% VNR0 = round(1e6*VNR0)/1e6

% fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180601/MATLAB/';
% fn = ['NR3-' erase(fn,'.txt') '-matlab' '.txt'];
% 
% fid = fopen([fp fn],'w');
% 
% fprintf(fid,'VNR\n\n');
% for ph = 1:3
%     fprintf(fid,[phstr{ph} '\n']);
%     for k1 = 1:nnode
%         fprintf(fid,'%1.6f %1.6f\n',real(VNR0(ph,k1)),imag(VNR0(ph,k1)));
%     end
% end
% 
% fprintf(fid,'\n\nINR\n\n');
% for ph = 1:3
%     fprintf(fid,[phstr{ph} '\n']);
%     for k1 = 1:nline
%         fprintf(fid,'%1.6f %1.6f\n',real(INR0(ph,k1)),imag(INR0(ph,k1)));
%     end
% end
% 
% fprintf(fid,'\n\nSTXNR\n\n');
% for ph = 1:3
%     fprintf(fid,[phstr{ph} '\n']);
%     for k1 = 1:nline
%         fprintf(fid,'%1.6f %1.6f\n',real(STXNR0(ph,k1)),imag(STXNR0(ph,k1)));
%     end
% end
% 
% fprintf(fid,'\n\nSRXNR\n\n');
% for ph = 1:3
%     fprintf(fid,[phstr{ph} '\n']);
%     for k1 = 1:nline
%         fprintf(fid,'%1.6f %1.6f\n',real(SRXNR0(ph,k1)),imag(SRXNR0(ph,k1)));
%     end
% end
% 
% fprintf(fid,'\n\niNR\n\n');
% for ph = 1:3
%     fprintf(fid,[phstr{ph} '\n']);
%     for k1 = 1:nnode
%         fprintf(fid,'%1.6f %1.6f\n',real(iNR0(ph,k1)),imag(iNR0(ph,k1)));
%     end
% end
% 
% fprintf(fid,'\n\nsNR\n\n');
% for ph = 1:3
%     fprintf(fid,[phstr{ph} '\n']);
%     for k1 = 1:nnode
%         fprintf(fid,'%1.6f %1.6f\n',real(sNR0(ph,k1)),imag(sNR0(ph,k1)));
%     end
% end
% 
% fclose(fid);