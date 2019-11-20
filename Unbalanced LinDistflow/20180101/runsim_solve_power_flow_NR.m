% Michael Sankur - msankur@berkeley.edu
% 2018.01.01

clc, clear, close all

% if exist('FBS_3phase_fun_20160603','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
%     path(path,'C:\Users\Michael\Desktop\reconfiguration\20160603');
% end

% path(path,genpath('C:\Users\Michael\Dropbox\Unbalanced LinDistflow\20180101\'));
path(path,genpath('/home/michael/Dropbox/Unbalanced LinDistflow/20180101/'));
% path(path,'C:\Users\Michael\Desktop\reconfiguration\Mesh\SDP');
% path(path,'C:\Users\Michael\Desktop\reconfiguration\Mesh\TX');

% if exist('mxlpsolve','file') == 0
%     path(path,'C:\Users\Michael\Desktop\mxlp');
% end
% 
% if exist('cvx_begin','file') == 0
%     cd C:\Users\Michael\Desktop\cvx-w64\cvx
%     cvx_setup
% end

%% Load feeder

feedername = 'ieee_13node_balance';
% feedername = 'ieee_37node_mesh_open';
% feedername = 'mesh_09node_balance';
% feedername = '05node_fullphase_radial';

% home = pwd; %returns a string for the path of the base folder in which the .m file is housed.
% if ispc %This should identify the appropriate pathseparator for the operating system so that it will run on a windows, mac, or linux machine.
%     pathseparator = '\';
% else
%     pathseparator = '/';
% end
% genpathstr=strcat(home); %path string to the GEN folder
% feederpathstr=strcat(home,pathseparator, 'Networks', pathseparator); %path string to the Networks folder
% path(path,genpath(genpathstr));

fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180101/NETWORKS/';
fn = [feedername '.txt'];

[network1] = network_mapper_function_20180101(feedername, fp, fn);

%% Network paramaters

nnode = network1.nodes.nnode;
nline = network1.lines.nline;

%% Load parameters

network1.loads.aPQ = 0.85*ones(3,nnode).*(network1.loads.spu ~= 0);
network1.loads.aI = zeros(3,nnode).*(network1.loads.spu ~= 0);
network1.loads.aZ = 0.15*ones(3,nnode).*(network1.loads.spu ~= 0);

%% Capacitor parameters

% network1.caps.cappu = 1*network1.caps.cappu;

%% Controller parameters

% network1.cons.wmaxpu = 1*network1.cons.wmaxpu;

%% IEEE 13 node params


% network1.loads.spu = 1.125*network1.loads.spu;

% network1.caps.cappu = 1*network1.caps.cappu;

% network1.cons.wmaxpu = 0.25*network1.cons.wmaxpu;

rho = 0.5;

%% Simulation parameters

slacknode = 1;
Vslack = 1*[1;
    1*exp(1j*-120*pi/180);
    1*exp(1j*120*pi/180)];

%%

network1.cons.wpu = zeros(3,nnode);

[VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, iter] = NR3(network1,[],[],slacknode,Vslack)

%%

phstr = {'a','b','c'};

% VNR0 = round(1e6*VNR0)/1e6

fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180101/';
fn = ['NR3' fn '-matlab' '.txt'];

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