% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

if exist('FBS_3phase_fun_20160603','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
    path(path,'C:\Users\Michael\Desktop\reconfiguration\20160603');
end

if exist('mxlpsolve','file') == 0
    path(path,'C:\Users\Michael\Desktop\mxlp');
end

if exist('cvx_begin','file') == 0
    cd C:\Users\Michael\Desktop\cvx-w64\cvx
    cvx_setup
end

%% Load feeder

feedername = 'ieee_13node_balance';
% feedername = '03node_fullphase';
feedername = '04node_fullphase';
feedername = '04node_multiphase';
feedername = '03node_fullphase';
feedername = '04node_singlephase';
feedername = '05node_singlephase_radial';
% feedername = '05node_singlephase_mesh';
% feedername = 'ieee_37node_singlephase';

fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fn = [feedername '.txt'];

[feeder, nodes, lines, configs, loads, caps, controllers] = feeder_mapper_TX_function_20170215(feedername, fp, fn);

%% Feeder paramaters

radflag1 = 1;
for k1 = 1:size(nodes.FM,1)
    if sum(nodes.FM(k1,:) == -1) > 1
        radflag1 = 0;
    end
end
radflag1

nnode = nodes.nnode;
nline = lines.nline;

%% Nodes parameters

%% Line parameters

%% Load parameters

% loads.aPQ = 0.85*ones(3,nnode).*nodes.PH;
% loads.aI = 0.0*ones(3,nnode).*nodes.PH;
% loads.aZ = 0.15*ones(3,nnode).*nodes.PH;

loads.spu = 1.5*loads.spu;

%% Capacitor parameters

caps.cappu = 1.0*caps.cappu;

%% Controller parameters

controllers.wmaxpu = 0.5*controllers.wmaxpu;

%% Simulation parameters

Vnom = 1;
sim.Vnom = Vnom;

slacknode = 1;
sim.slacknode = 1;

sim.rho = 0.1;

%%

[sdpsol,SDPMAT] = Solver_SDP_minP0_singlephase(feeder, nodes, lines, configs, loads, caps, controllers, sim)

[V, D] = eig(sdpsol.Xsdp)

%%

Vnom = [1;
    1*exp(j*240*pi/180);
    1*exp(j*120*pi/180)];

controllers.wpu = zeros(3,nnode);

[VNR0, INR0, STXNR0, SRXNR0] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,[],[],slacknode,Vnom);

controllers.wpu = sdpsol.wsdp;
controllers.wpu(:,1) = 0;

[VNR1, INR1, STXNR1, SRXNR1] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,[],[],slacknode,Vnom);

%%

