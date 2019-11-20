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

loads.origspu = loads.spu;

ma = 0.5:0.025:2.0;

EERROR = [];
DERROR = [];
VERROR = [];
PERROR = [];
QERROR = [];
SERROR = [];

for k1 = 1:length(ma)
    
    loads.spu = ma(k1)*loads.origspu;
    
    FBS0 = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

    FBS1 = FBS_TX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

    FBS2 = FBS_TX_A2(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

    FBS3 = FBS_TX_A3(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

    FBS4 = FBS_TX_A4(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

    FBS5 = FBS_TX_A5(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

    FBS6 = FBS_TX_A6(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    
    EERROR(:,k1) = [max(max(abs(abs(FBS0.V) - abs(FBS1.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS2.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS3.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS4.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS5.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS6.V))))];

    DERROR(:,k1) = [max(max(abs(angle(FBS0.V) - angle(FBS1.V))));
        max(max(abs(angle(FBS0.V) - angle(FBS2.V))));
        max(max(abs(angle(FBS0.V) - angle(FBS3.V))));
        max(max(abs(angle(FBS0.V) - angle(FBS4.V))));
        max(max(abs(angle(FBS0.V) - angle(FBS5.V))));
        max(max(abs(angle(FBS0.V) - angle(FBS6.V))))];

    VERROR(:,k1) = [max(max(abs(FBS0.V - FBS1.V)));
        max(max(abs(FBS0.V - FBS2.V)));
        max(max(abs(FBS0.V - FBS3.V)));
        max(max(abs(FBS0.V - FBS4.V)));
        max(max(abs(FBS0.V - FBS5.V)));
        max(max(abs(FBS0.V - FBS6.V)))];

    PERROR(:,k1) = [max(max(abs(real(FBS0.Srx - FBS1.Srx))));
        max(max(abs(real(FBS0.Srx - FBS2.Srx))));
        max(max(abs(real(FBS0.Srx - FBS3.Srx))));
        max(max(abs(real(FBS0.Srx - FBS4.Srx))));
        max(max(abs(real(FBS0.Srx - FBS5.Srx))));
        max(max(abs(real(FBS0.Srx - FBS6.Srx))))];

    QERROR(:,k1) = [max(max(abs(imag(FBS0.Srx - FBS1.Srx))));
        max(max(abs(imag(FBS0.Srx - FBS2.Srx))));
        max(max(abs(imag(FBS0.Srx - FBS3.Srx))));
        max(max(abs(imag(FBS0.Srx - FBS4.Srx))));
        max(max(abs(imag(FBS0.Srx - FBS5.Srx))));
        max(max(abs(imag(FBS0.Srx - FBS6.Srx))))];

    SERROR(:,k1) = [max(max(abs(FBS0.Srx - FBS1.Srx)));
        max(max(abs(FBS0.Srx - FBS2.Srx)));
        max(max(abs(FBS0.Srx - FBS3.Srx)));
        max(max(abs(FBS0.Srx - FBS4.Srx)));
        max(max(abs(FBS0.Srx - FBS5.Srx)));
        max(max(abs(FBS0.Srx - FBS6.Srx)))];
    
    
end

%%

lss = {'r.','g.','b.','m.','k.','c.'};

figure, box on, hold on
plot(ma,EERROR(1,:),lss{1},'LineWidth',2)
plot(ma,EERROR(2,:),lss{2},'LineWidth',2)
plot(ma,EERROR(3,:),lss{3},'LineWidth',2)
plot(ma,EERROR(4,:),lss{4},'LineWidth',2)
plot(ma,EERROR(5,:),lss{5},'LineWidth',2)
plot(ma,EERROR(6,:),lss{6},'LineWidth',2)

figure, box on, hold on
plot(ma,DERROR(1,:),lss{1},'LineWidth',2)
plot(ma,DERROR(2,:),lss{2},'LineWidth',2)
plot(ma,DERROR(3,:),lss{3},'LineWidth',2)
plot(ma,DERROR(4,:),lss{4},'LineWidth',2)
plot(ma,DERROR(5,:),lss{5},'LineWidth',2)
plot(ma,DERROR(6,:),lss{6},'LineWidth',2)

figure, box on, hold on
plot(ma,VERROR(1,:),lss{1},'LineWidth',2)
plot(ma,VERROR(2,:),lss{2},'LineWidth',2)
plot(ma,VERROR(3,:),lss{3},'LineWidth',2)
plot(ma,VERROR(4,:),lss{4},'LineWidth',2)
plot(ma,VERROR(5,:),lss{5},'LineWidth',2)
plot(ma,VERROR(6,:),lss{6},'LineWidth',2)

figure, box on, hold on
plot(ma,PERROR(1,:),lss{1},'LineWidth',2)
plot(ma,PERROR(2,:),lss{2},'LineWidth',2)
plot(ma,PERROR(3,:),lss{3},'LineWidth',2)
plot(ma,PERROR(4,:),lss{4},'LineWidth',2)
plot(ma,PERROR(5,:),lss{5},'LineWidth',2)
plot(ma,PERROR(6,:),lss{6},'LineWidth',2)

figure, box on, hold on
plot(ma,QERROR(1,:),lss{1},'LineWidth',2)
plot(ma,QERROR(2,:),lss{2},'LineWidth',2)
plot(ma,QERROR(3,:),lss{3},'LineWidth',2)
plot(ma,QERROR(4,:),lss{4},'LineWidth',2)
plot(ma,QERROR(5,:),lss{5},'LineWidth',2)
plot(ma,QERROR(6,:),lss{6},'LineWidth',2)

figure, box on, hold on
plot(ma,SERROR(1,:),lss{1},'LineWidth',2)
plot(ma,SERROR(2,:),lss{2},'LineWidth',2)
plot(ma,SERROR(3,:),lss{3},'LineWidth',2)
plot(ma,SERROR(4,:),lss{4},'LineWidth',2)
plot(ma,SERROR(5,:),lss{5},'LineWidth',2)
plot(ma,SERROR(6,:),lss{6},'LineWidth',2)





%%

% FBS0 = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
% FBS1 = FBS_TX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
% FBS2 = FBS_TX_A2(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
% FBS3 = FBS_TX_A3(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
% FBS4 = FBS_TX_A4(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
% FBS5 = FBS_TX_A5(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
% FBS6 = FBS_TX_A6(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
% EERROR = [max(max(abs(abs(FBS0.V) - abs(FBS1.V))));
%     max(max(abs(abs(FBS0.V) - abs(FBS2.V))));
%     max(max(abs(abs(FBS0.V) - abs(FBS3.V))));
%     max(max(abs(abs(FBS0.V) - abs(FBS4.V))));
%     max(max(abs(abs(FBS0.V) - abs(FBS5.V))));
%     max(max(abs(abs(FBS0.V) - abs(FBS6.V))))];
% 
% DERROR = [max(max(abs(angle(FBS0.V) - angle(FBS1.V))));
%     max(max(abs(angle(FBS0.V) - angle(FBS2.V))));
%     max(max(abs(angle(FBS0.V) - angle(FBS3.V))));
%     max(max(abs(angle(FBS0.V) - angle(FBS4.V))));
%     max(max(abs(angle(FBS0.V) - angle(FBS5.V))));
%     max(max(abs(angle(FBS0.V) - angle(FBS6.V))))];
% 
% VERROR = [max(max(abs(FBS0.V - FBS1.V)));
%     max(max(abs(FBS0.V - FBS2.V)));
%     max(max(abs(FBS0.V - FBS3.V)));
%     max(max(abs(FBS0.V - FBS4.V)));
%     max(max(abs(FBS0.V - FBS5.V)));
%     max(max(abs(FBS0.V - FBS6.V)))];
% 
% PERROR = [max(max(abs(real(FBS0.Srx - FBS1.Srx))));
%     max(max(abs(real(FBS0.Srx - FBS2.Srx))));
%     max(max(abs(real(FBS0.Srx - FBS3.Srx))));
%     max(max(abs(real(FBS0.Srx - FBS4.Srx))));
%     max(max(abs(real(FBS0.Srx - FBS5.Srx))));
%     max(max(abs(real(FBS0.Srx - FBS6.Srx))))];
% 
% QERROR = [max(max(abs(imag(FBS0.Srx - FBS1.Srx))));
%     max(max(abs(imag(FBS0.Srx - FBS2.Srx))));
%     max(max(abs(imag(FBS0.Srx - FBS3.Srx))));
%     max(max(abs(imag(FBS0.Srx - FBS4.Srx))));
%     max(max(abs(imag(FBS0.Srx - FBS5.Srx))));
%     max(max(abs(imag(FBS0.Srx - FBS6.Srx))))];
% 
% SERROR = [max(max(abs(FBS0.Srx - FBS1.Srx)));
%     max(max(abs(FBS0.Srx - FBS2.Srx)));
%     max(max(abs(FBS0.Srx - FBS3.Srx)));
%     max(max(abs(FBS0.Srx - FBS4.Srx)));
%     max(max(abs(FBS0.Srx - FBS5.Srx)));
%     max(max(abs(FBS0.Srx - FBS6.Srx)))];

