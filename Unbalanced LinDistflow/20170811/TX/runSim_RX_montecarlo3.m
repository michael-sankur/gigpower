% Michael Sankur - msankur@lbl.gov
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

if exist('FBS','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
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

% name = '6node_fullphase_test';
name = 'ieee_13node_new';

fp = ['C:\Users\Michael\Desktop\reconfiguration\Feeders\'];
fn = [name '.txt'];

[feeder, nodes, lines, configs, loads, caps, controllers] = feeder_mapper_TX_function_20170215(name, fp, fn);


%% Feeder paramaters

nnode = nodes.nnode;
nline = lines.nline;

%% Load parameters

loads.aPQ = 0.85*ones(3,nnode).*nodes.PH;
loads.aI = zeros(3,nnode);
loads.aZ = 0.15*ones(3,nnode).*nodes.PH;

loads.spu = 1.0*loads.spu;

%% Controller parameters

controllers.wmaxpu = 1*controllers.wmaxpu;

controllers.wpu = zeros(3,nnode);

%% FBS New

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim2.Vnom = Vnom;

sim2.Vfbs = [];
sim2.Lfbs = [];
sim2.Hfbs = [];

sim.rho = 0.1;

%% Monte Carlo

lines.origFZpu = lines.FZpu;
loads.origspu = loads.spu;

ma1 = 0.5:0.025:2;
ma2 = 0.75:0.025:1.25;

EERROR = [];
DERROR = [];
VERROR = [];
PERROR = [];
QERROR = [];
SERROR = [];

for k1 = 1:length(ma1)
    for k2 = 1:length(ma2)
    
        loads.spu = ma1(k1)*loads.origspu;
        lines.FZpu = ma2(k2)*lines.origFZpu;

        FBS0 = FBS_RX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

        FBS1 = FBS_RX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

%         FBS2 = FBS_RX_A2(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
%         FBS3 = FBS_RX_A3(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
%         FBS4 = FBS_RX_A4(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
%         FBS5 = FBS_RX_A5(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
%         FBS6 = FBS_RX_A6(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

        EERROR(k1,k2) = max(max(abs(abs(FBS0.V) - abs(FBS1.V))));
        DERROR(k1,k2) = max(max(abs(angle(FBS0.V) - angle(FBS1.V))));
        VERROR(k1,k2) = max(max(abs(FBS0.V - FBS1.V)));
        
        PERROR(k1,k2) = max(max(abs(real(FBS0.Srx - FBS1.Srx))));
        QERROR(k1,k2) = max(max(abs(imag(FBS0.Srx - FBS1.Srx))));
        SERROR(k1,k2) = max(max(abs(FBS0.Srx - FBS1.Srx)));

    end
end

%%

lss = {'r.','g.','b.','m.','k.','c.'};

figure, box on, hold on
% surf(ma1,ma2,EERROR','LineStyle','none')
surf(ma1,ma2,DERROR','LineStyle','none')
