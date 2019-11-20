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

name = '9node_fullphase_test';
% name = 'ieee_13node_new';


fp = ['C:\Users\Michael\Desktop\reconfiguration\Feeders\'];
fn = [name '.txt'];

[feeder, nodes, lines, configs, loads, caps, controllers] = feeder_mapper_TX_function_20170215(name, fp, fn);


%% Feeder paramaters

nnode = nodes.nnode;
nline = lines.nline;

%% Line parameters and adding extra line

% lines.TXnum = [0 lines.TXnum];
% lines.TXname = ['SUB' lines.TXname];
% 
% lines.RXnum = [1 lines.RXnum];
% lines.RXname = ['A1' lines.RXname];
% 
% lines.PH = [ones(3,1) lines.PH];
% 
% lines.config = ['12' lines.config];
% 
% lines.length = [3.048 lines.length];
% 
% lines.nline = lines.nline + 1;
% 
% lines.FZ(:,:,2:end+1) = lines.FZ;
% lines.FZ(:,:,1) = zeros(3,3);
% 
% lines.FZpu(:,:,2:end+1) = lines.FZpu;
% lines.FZpu(:,:,1) = zeros(3,3);
% 
% nline = lines.nline;

%% Load parameters

loads.spu = 2*loads.spu;

%% Controller parameters

controllers.wmaxpu = 0*controllers.wmaxpu

%% Simulation parameters

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim.Vnom = Vnom;

sim.Vfbs = [];
sim.Lfbs = [];
sim.Hfbs = [];

sim.rho = 0.1;

%%

OPTSOL = Solver_LinDist3Flow_TX_minP0_20160603_mag(feeder, nodes, lines, configs, loads, caps, controllers, sim)

nvar = OPTSOL.nvar;
Xa = OPTSOL.X(1:nvar);
Xb = OPTSOL.X(nvar+1:2*nvar);
Xc = OPTSOL.X(2*nvar+1:3*nvar);

XX = [Xa Xb Xc]';

Yopt1 = XX(:,1:nnode);
XX(:,1:nnode) = [];

Yopt1
Vopt1 = sqrt(Yopt1)

Popt1 = XX(:,1:nline);
XX(:,1:nline) = [];

Qopt1 = XX(:,1:nline);
XX(:,1:nline) = [];

Sopt1 = Popt1 + 1j*Qopt1

uopt1 = XX(:,1:nnode);
XX(:,1:nnode) = [];

vopt1 = XX(:,1:nnode);
XX(:,1:nnode) = [];

wopt1 = uopt1 + 1j*vopt1

dem1 = loads.spu.*(loads.aPQ + loads.aZ.*Yopt1).*nodes.PH

for k1 = 2:nnode
    sopt1(:,k1) = sum(Sopt1(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
        sum(Sopt1(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2);
end
sopt1

%%

OPTSOL = Solver_LinDist3Flow_TX_minP0_20160603_magangle(feeder, nodes, lines, configs, loads, caps, controllers, sim)

nvar = OPTSOL.nvar;
Xa = OPTSOL.X(1:nvar);
Xb = OPTSOL.X(nvar+1:2*nvar);
Xc = OPTSOL.X(2*nvar+1:3*nvar);

XX = [Xa Xb Xc]';

Yopt2 = XX(:,1:nnode);
XX(:,1:nnode) = [];

Yopt2
Vopt2 = sqrt(Yopt2)

Dopt2 = XX(:,1:nnode);
XX(:,1:nnode) = [];

Dopt2

Popt2 = XX(:,1:nline);
XX(:,1:nline) = [];

Qopt2 = XX(:,1:nline);
XX(:,1:nline) = [];

Sopt2 = Popt2 + 1j*Qopt2

uopt2 = XX(:,1:nnode);
XX(:,1:nnode) = [];

vopt2 = XX(:,1:nnode);
XX(:,1:nnode) = [];

wopt2 = uopt2 + 1j*vopt2

dem2 = loads.spu.*(loads.aPQ + loads.aZ.*Yopt2).*nodes.PH

for k1 = 2:nnode
    sopt2(:,k1) = sum(Sopt2(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
        sum(Sopt2(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2);
end
sopt2

%%

Yopt1 - Yopt2
Vopt1 - Vopt2
Sopt1 - Sopt2
wopt1 - wopt2
dem1 - dem2
sopt1 - sopt2

%%

OPTSOL = Solver_LinDist3Flow_TX_minSmn_20160603_mag(feeder, nodes, lines, configs, loads, caps, controllers, sim)

nvar = OPTSOL.nvar;
Xa = OPTSOL.X(1:nvar);
Xb = OPTSOL.X(nvar+1:2*nvar);
Xc = OPTSOL.X(2*nvar+1:3*nvar);

XX = [Xa Xb Xc]';

Yopt3 = XX(:,1:nnode);
XX(:,1:nnode) = [];

Yopt3
Vopt3 = sqrt(Yopt3)

Popt3 = XX(:,1:nline);
XX(:,1:nline) = [];

Qopt3 = XX(:,1:nline);
XX(:,1:nline) = [];

Sopt3 = Popt3 + 1j*Qopt3

uopt3 = XX(:,1:nnode);
XX(:,1:nnode) = [];

vopt3 = XX(:,1:nnode);
XX(:,1:nnode) = [];

wopt3 = uopt3 + 1j*vopt3

dem3 = loads.spu.*(loads.aPQ + loads.aZ.*Yopt3).*nodes.PH

for k1 = 2:nnode
    sopt3(:,k1) = sum(Sopt3(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
        sum(Sopt3(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2);
end
sopt3

%%

OPTSOL = Solver_LinDist3Flow_TX_minSmn_20160603_magangle(feeder, nodes, lines, configs, loads, caps, controllers, sim)

nvar = OPTSOL.nvar;
Xa = OPTSOL.X(1:nvar);
Xb = OPTSOL.X(nvar+1:2*nvar);
Xc = OPTSOL.X(2*nvar+1:3*nvar);

XX = [Xa Xb Xc]';

Yopt4 = XX(:,1:nnode);
XX(:,1:nnode) = [];

Yopt4
Vopt4 = sqrt(Yopt4)

Dopt4 = XX(:,1:nnode);
XX(:,1:nnode) = [];

Dopt4

Popt4 = XX(:,1:nline);
XX(:,1:nline) = [];

Qopt4 = XX(:,1:nline);
XX(:,1:nline) = [];

Sopt4 = Popt4 + 1j*Qopt4

uopt4 = XX(:,1:nnode);
XX(:,1:nnode) = [];

vopt4 = XX(:,1:nnode);
XX(:,1:nnode) = [];

wopt4 = uopt4 + 1j*vopt4

dem4 = loads.spu.*(loads.aPQ + loads.aZ.*Yopt4).*nodes.PH

for k1 = 2:nnode
    sopt4(:,k1) = sum(Sopt4(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
        sum(Sopt4(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2);
end
sopt4

%%

Yopt3 - Yopt4
Vopt3 - Vopt4
Sopt3 - Sopt4
wopt3 - wopt4
dem3 - dem4
sopt3 - sopt4
