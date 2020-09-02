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
% feedername = '03node_fullphase';
% feedername = '04node_singlephase';

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

%% Line parameters

lines.FYpu(:,:,1) = 1*lines.FYpu(:,:,1)

%% Load parameters

loads.aPQ = 1.0*ones(3,nnode).*nodes.PH;
loads.aI = 0.0*ones(3,nnode).*nodes.PH;
loads.aZ = 0.0*ones(3,nnode).*nodes.PH;

loads.spu = 1.5*loads.spu;

%% Capacitor parameters

caps.cappu = 1.0*caps.cappu;

%% Controller parameters

controllers.wmaxpu = 1.0*controllers.wmaxpu;

%% Simulation parameters

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim.Vnom = Vnom;

sim.Vfbs = [];
sim.Lfbs = [];
sim.Hfbs = [];

slacknode = 1;
sim.slacknode = 1;

sim.rho = 0.1;

%%

[sdpsol,SDPMAT] = Solver_SDP_TX_minP0(feeder, nodes, lines, configs, loads, caps, controllers, sim)

%%

controllers.wpu = zeros(3,nnode);

[VNR0, INR0, STXNR0, SRXNR0] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,[],[],slacknode,Vnom);

% sdpsol.wsdp(:,1) = 0;
controllers.wpu = sdpsol.wsdp;

[VNR1, INR1, STXNR1, SRXNR1] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,[],[],slacknode,Vnom);

%%

BCIB = NaN*ones(1,nnode);
CCIB = NaN*ones(1,nnode);

for k1 = 2:nnode
    if nodes.PH(1,k1) == 1 && nodes.PH(2,k1) == 1 && nodes.PH(3,k1) == 0
        BCIB(k1) = abs(abs(VNR0(1,k1)) - abs(VNR0(2,k1)));
        CCIB(k1) = abs(abs(VNR1(1,k1)) - abs(VNR1(2,k1)));
    elseif nodes.PH(1,k1) == 1 && nodes.PH(2,k1) == 0 && nodes.PH(3,k1) == 1
        BCIB(k1) = abs(abs(VNR0(1,k1)) - abs(VNR0(3,k1)));
        CCIB(k1) = abs(abs(VNR1(1,k1)) - abs(VNR1(3,k1)));
    elseif nodes.PH(1,k1) == 0 && nodes.PH(2,k1) == 1 && nodes.PH(3,k1) == 1
        BCIB(k1) = abs(abs(VNR0(2,k1)) - abs(VNR0(3,k1)));
        CCIB(k1) = abs(abs(VNR1(2,k1)) - abs(VNR1(3,k1)));
    elseif nodes.PH(1,k1) == 1 && nodes.PH(2,k1) == 1 && nodes.PH(3,k1) == 1
        BCIB(k1) = abs(abs(VNR0(1,k1)) - abs(VNR0(2,k1))) ...
            + abs(abs(VNR0(1,k1)) - abs(VNR0(3,k1))) ...
            + abs(abs(VNR0(2,k1)) - abs(VNR0(3,k1)));
        CCIB(k1) = abs(abs(VNR1(1,k1)) - abs(VNR1(2,k1))) ...
            + abs(abs(VNR1(1,k1)) - abs(VNR1(3,k1))) ...
            + abs(abs(VNR1(2,k1)) - abs(VNR1(3,k1)));        
    end    
    
end

sum(BCIB(isnan(BCIB) == 0))
sum(CCIB(isnan(CCIB) == 0))


%%

% nodes.nodelist{1} = 'TXL';

close all

VNR0PLOT = VNR0; VNR0PLOT(nodes.PH == 0) = NaN;
VNR1PLOT = VNR1; VNR1PLOT(nodes.PH == 0) = NaN;

figure, box on, hold on
plot(1:nnode,abs(VNR0PLOT(1,:)),'r+','MarkerSize',12.5,'LineWidth',2)
plot(1:nnode,abs(VNR0PLOT(2,:)),'gx','MarkerSize',12.5,'LineWidth',2)
plot(1:nnode,abs(VNR0PLOT(3,:)),'b.','MarkerSize',25,'LineWidth',2)
plot(0:nnode+1,0.95*ones(size(0:nnode+1)),'k--','LineWidth',2)
set(gca,'XTick',1:nnode,'XTickLabel',nodes.nodelist,'YTick',0.94:0.01:1.01, ...
    'FontWeight','bold','FontSize',15)
legend({'A','B','C'},'FontWeight','bold','FontSize',12)
title('Base Case - No Control','FontWeight','bold','FontSize',15)
xlabel('Node','FontWeight','bold','FontSize',12)
ylabel('Voltage Magnitude [p.u.]','FontWeight','bold','FontSize',12)
axis([0.5 nnode+0.5 0.94 1.01])
pbaspect([1 0.5 1])

print('-f1','-depsc','C:\Users\Michael\Desktop\temp\eps\balance13nodebase.eps')


figure, box on, hold on
plot(1:nnode,abs(VNR1PLOT(1,:)),'r+','MarkerSize',12.5,'LineWidth',2)
plot(1:nnode,abs(VNR1PLOT(2,:)),'gx','MarkerSize',12.5,'LineWidth',2)
plot(1:nnode,abs(VNR1PLOT(3,:)),'b.','MarkerSize',25,'LineWidth',2)
plot(0:nnode+1,0.95*ones(size(0:nnode+1)),'k--','LineWidth',2)
set(gca,'XTick',1:nnode,'XTickLabel',nodes.nodelist,'YTick',0.94:0.01:1.01, ...
    'FontWeight','bold','FontSize',15)
legend({'A','B','C'},'FontWeight','bold','FontSize',12)
title('Control Case','FontWeight','bold','FontSize',15)
xlabel('Node','FontWeight','bold','FontSize',12)
ylabel('Voltage Magnitude [p.u.]','FontWeight','bold','FontSize',12)
axis([0.5 nnode+0.5 0.94 1.01])
pbaspect([1 0.5 1])

print('-f2','-depsc','C:\Users\Michael\Desktop\temp\eps\balance13nodecontrol.eps')


figure, box on, hold on
plot(1:nnode,BCIB,'r+','MarkerSize',12.5,'LineWidth',2)
plot(1:nnode,CCIB,'gx','MarkerSize',12.5,'LineWidth',2)
set(gca,'XTick',1:nnode,'XTickLabel',nodes.nodelist, ...
    'FontWeight','bold','FontSize',15)
legend({'Base Case','Control Case'},'FontWeight','bold','FontSize',12,'location','northwest')
title('Voltage Imbalance','FontWeight','bold','FontSize',15)
xlabel('Node','FontWeight','bold','FontSize',12)
ylabel('Voltage Imbalance [p.u.]','FontWeight','bold','FontSize',12)
axis([0.5 nnode+0.5 0 0.1])
pbaspect([1 0.5 1])

print('-f3','-depsc','C:\Users\Michael\Desktop\temp\eps\balance13nodeimbalance.eps')




