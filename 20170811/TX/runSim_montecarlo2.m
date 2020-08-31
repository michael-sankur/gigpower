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
    
    FBSTX0 = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    FBSTX1 = FBS_TX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    FBSTX2 = FBS_TX_A2(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    FBSTX3 = FBS_TX_A3(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    FBSTX4 = FBS_TX_A4(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    FBSTX5 = FBS_TX_A5(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    FBSTX6 = FBS_TX_A6(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    
    FBSRX0 = FBS_RX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);   
    FBSRX1 = FBS_RX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);    
    FBSRX2 = FBS_RX_A2(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    FBSRX3 = FBS_RX_A3(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    FBSRX4 = FBS_RX_A4(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);    
    FBSRX5 = FBS_RX_A5(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);    
    FBSRX6 = FBS_RX_A6(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
    
    TXEERROR(:,k1) = [max(max(abs(abs(FBSTX0.V) - abs(FBSTX1.V))));
        max(max(abs(abs(FBSTX0.V) - abs(FBSTX2.V))));
        max(max(abs(abs(FBSTX0.V) - abs(FBSTX3.V))));
        max(max(abs(abs(FBSTX0.V) - abs(FBSTX4.V))));
        max(max(abs(abs(FBSTX0.V) - abs(FBSTX5.V))));
        max(max(abs(abs(FBSTX0.V) - abs(FBSTX6.V))))];
    
    RXEERROR(:,k1) = [max(max(abs(abs(FBSRX0.V) - abs(FBSRX1.V))));
        max(max(abs(abs(FBSRX0.V) - abs(FBSRX2.V))));
        max(max(abs(abs(FBSRX0.V) - abs(FBSRX3.V))));
        max(max(abs(abs(FBSRX0.V) - abs(FBSRX4.V))));
        max(max(abs(abs(FBSRX0.V) - abs(FBSRX5.V))));
        max(max(abs(abs(FBSRX0.V) - abs(FBSRX6.V))))];

    TXDERROR(:,k1) = [max(max(abs(angle(FBSTX0.V) - angle(FBSTX1.V))));
        max(max(abs(angle(FBSTX0.V) - angle(FBSTX2.V))));
        max(max(abs(angle(FBSTX0.V) - angle(FBSTX3.V))));
        max(max(abs(angle(FBSTX0.V) - angle(FBSTX4.V))));
        max(max(abs(angle(FBSTX0.V) - angle(FBSTX5.V))));
        max(max(abs(angle(FBSTX0.V) - angle(FBSTX6.V))))];
    
    RXDERROR(:,k1) = [max(max(abs(angle(FBSRX0.V) - angle(FBSRX1.V))));
        max(max(abs(angle(FBSRX0.V) - angle(FBSRX2.V))));
        max(max(abs(angle(FBSRX0.V) - angle(FBSRX3.V))));
        max(max(abs(angle(FBSRX0.V) - angle(FBSRX4.V))));
        max(max(abs(angle(FBSRX0.V) - angle(FBSRX5.V))));
        max(max(abs(angle(FBSRX0.V) - angle(FBSRX6.V))))];
% 
%     VERROR(:,k1) = [max(max(abs(FBSTX0.V - FBSTX1.V)));
%         max(max(abs(FBSTX0.V - FBSTX2.V)));
%         max(max(abs(FBSTX0.V - FBSTX3.V)));
%         max(max(abs(FBSTX0.V - FBSTX4.V)));
%         max(max(abs(FBSTX0.V - FBSTX5.V)));
%         max(max(abs(FBSTX0.V - FBSTX6.V)))];
% 
%     PERROR(:,k1) = [max(max(abs(real(FBSTX0.Srx - FBSTX1.Srx))));
%         max(max(abs(real(FBSTX0.Srx - FBSTX2.Srx))));
%         max(max(abs(real(FBSTX0.Srx - FBSTX3.Srx))));
%         max(max(abs(real(FBSTX0.Srx - FBSTX4.Srx))));
%         max(max(abs(real(FBSTX0.Srx - FBSTX5.Srx))));
%         max(max(abs(real(FBSTX0.Srx - FBSTX6.Srx))))];
% 
%     QERROR(:,k1) = [max(max(abs(imag(FBSTX0.Srx - FBSTX1.Srx))));
%         max(max(abs(imag(FBSTX0.Srx - FBSTX2.Srx))));
%         max(max(abs(imag(FBSTX0.Srx - FBSTX3.Srx))));
%         max(max(abs(imag(FBSTX0.Srx - FBSTX4.Srx))));
%         max(max(abs(imag(FBSTX0.Srx - FBSTX5.Srx))));
%         max(max(abs(imag(FBSTX0.Srx - FBSTX6.Srx))))];
% 
    TXSERROR(:,k1) = [max(max(abs(FBSTX0.Srx - FBSTX1.Srx)));
        max(max(abs(FBSTX0.Srx - FBSTX2.Srx)));
        max(max(abs(FBSTX0.Srx - FBSTX3.Srx)));
        max(max(abs(FBSTX0.Srx - FBSTX4.Srx)));
        max(max(abs(FBSTX0.Srx - FBSTX5.Srx)));
        max(max(abs(FBSTX0.Srx - FBSTX6.Srx)))];
    
    RXSERROR(:,k1) = [max(max(abs(FBSRX0.Srx - FBSRX1.Srx)));
        max(max(abs(FBSRX0.Srx - FBSRX2.Srx)));
        max(max(abs(FBSRX0.Srx - FBSRX3.Srx)));
        max(max(abs(FBSRX0.Srx - FBSRX4.Srx)));
        max(max(abs(FBSRX0.Srx - FBSRX5.Srx)));
        max(max(abs(FBSRX0.Srx - FBSRX6.Srx)))];
    
    
end

%%

txlss = {'r.','g.','b.','m.','k.','c.'};
rxlss = {'ro','go','bo','mo','ko','co'};

figure, box on, hold on
plot(ma,TXEERROR(1,:),txlss{1},'LineWidth',2)
plot(ma,TXEERROR(2,:),txlss{2},'LineWidth',2)
plot(ma,TXEERROR(3,:),txlss{3},'LineWidth',2)
plot(ma,TXEERROR(4,:),txlss{4},'LineWidth',2)
plot(ma,TXEERROR(5,:),txlss{5},'LineWidth',2)
plot(ma,TXEERROR(6,:),txlss{6},'LineWidth',2)
plot(ma,RXEERROR(1,:),rxlss{1},'LineWidth',2)
plot(ma,RXEERROR(2,:),rxlss{2},'LineWidth',2)
plot(ma,RXEERROR(3,:),rxlss{3},'LineWidth',2)
plot(ma,RXEERROR(4,:),rxlss{4},'LineWidth',2)
plot(ma,RXEERROR(5,:),rxlss{5},'LineWidth',2)
plot(ma,RXEERROR(6,:),rxlss{6},'LineWidth',2)
set(gca,'Fontsize',15,'FontWeight','bold')
title('Voltage Magnitude Error','Fontsize',15,'FontWeight','bold')
xlabel('IEEE 13 node load multiplier','Fontsize',15,'FontWeight','bold')
ylabel('Voltage Error [p.u.]','Fontsize',15,'FontWeight','bold')


figure, box on, hold on
plot(ma,TXDERROR(1,:),txlss{1},'LineWidth',2)
plot(ma,TXDERROR(2,:),txlss{2},'LineWidth',2)
plot(ma,TXDERROR(3,:),txlss{3},'LineWidth',2)
plot(ma,TXDERROR(4,:),txlss{4},'LineWidth',2)
plot(ma,TXDERROR(5,:),txlss{5},'LineWidth',2)
plot(ma,TXDERROR(6,:),txlss{6},'LineWidth',2)
plot(ma,RXDERROR(1,:),rxlss{1},'LineWidth',2)
plot(ma,RXDERROR(2,:),rxlss{2},'LineWidth',2)
plot(ma,RXDERROR(3,:),rxlss{3},'LineWidth',2)
plot(ma,RXDERROR(4,:),rxlss{4},'LineWidth',2)
plot(ma,RXDERROR(5,:),rxlss{5},'LineWidth',2)
plot(ma,RXDERROR(6,:),rxlss{6},'LineWidth',2)

% figure, box on, hold on
% plot(ma,DERROR(1,:),lss{1},'LineWidth',2)
% plot(ma,DERROR(2,:),lss{2},'LineWidth',2)
% plot(ma,DERROR(3,:),lss{3},'LineWidth',2)
% plot(ma,DERROR(4,:),lss{4},'LineWidth',2)
% plot(ma,DERROR(5,:),lss{5},'LineWidth',2)
% plot(ma,DERROR(6,:),lss{6},'LineWidth',2)
% 
% figure, box on, hold on
% plot(ma,VERROR(1,:),lss{1},'LineWidth',2)
% plot(ma,VERROR(2,:),lss{2},'LineWidth',2)
% plot(ma,VERROR(3,:),lss{3},'LineWidth',2)
% plot(ma,VERROR(4,:),lss{4},'LineWidth',2)
% plot(ma,VERROR(5,:),lss{5},'LineWidth',2)
% plot(ma,VERROR(6,:),lss{6},'LineWidth',2)
% 
% figure, box on, hold on
% plot(ma,PERROR(1,:),lss{1},'LineWidth',2)
% plot(ma,PERROR(2,:),lss{2},'LineWidth',2)
% plot(ma,PERROR(3,:),lss{3},'LineWidth',2)
% plot(ma,PERROR(4,:),lss{4},'LineWidth',2)
% plot(ma,PERROR(5,:),lss{5},'LineWidth',2)
% plot(ma,PERROR(6,:),lss{6},'LineWidth',2)
% 
% figure, box on, hold on
% plot(ma,QERROR(1,:),lss{1},'LineWidth',2)
% plot(ma,QERROR(2,:),lss{2},'LineWidth',2)
% plot(ma,QERROR(3,:),lss{3},'LineWidth',2)
% plot(ma,QERROR(4,:),lss{4},'LineWidth',2)
% plot(ma,QERROR(5,:),lss{5},'LineWidth',2)
% plot(ma,QERROR(6,:),lss{6},'LineWidth',2)

figure, box on, hold on
plot(ma,TXSERROR(1,:),txlss{1},'LineWidth',2)
plot(ma,TXSERROR(2,:),txlss{2},'LineWidth',2)
plot(ma,TXSERROR(3,:),txlss{3},'LineWidth',2)
plot(ma,TXSERROR(4,:),txlss{4},'LineWidth',2)
plot(ma,TXSERROR(5,:),txlss{5},'LineWidth',2)
plot(ma,TXSERROR(6,:),txlss{6},'LineWidth',2)
plot(ma,RXSERROR(1,:),rxlss{1},'LineWidth',2)
plot(ma,RXSERROR(2,:),rxlss{2},'LineWidth',2)
plot(ma,RXSERROR(3,:),rxlss{3},'LineWidth',2)
plot(ma,RXSERROR(4,:),rxlss{4},'LineWidth',2)
plot(ma,RXSERROR(5,:),rxlss{5},'LineWidth',2)
plot(ma,RXSERROR(6,:),rxlss{6},'LineWidth',2)





