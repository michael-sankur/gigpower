% Michael Sankur - msankur@lbl.gov
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

if exist('FBSTX','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
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

%% FBSTX New

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim2.Vnom = Vnom;

sim2.Vfbs = [];
sim2.Lfbs = [];
sim2.Hfbs = [];

sim.rho = 0.1;

%%

FBSTX0 = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
FBSRX0 = FBS_RX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBSTX1 = FBS_TX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
FBSRX1 = FBS_RX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBSTX2 = FBS_TX_A2(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
FBSRX2 = FBS_RX_A2(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBSTX3 = FBS_TX_A3(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
FBSRX3 = FBS_RX_A3(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBSTX4 = FBS_TX_A4(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
FBSRX4 = FBS_RX_A4(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBSTX5 = FBS_TX_A5(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
FBSRX5 = FBS_RX_A5(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBSTX6 = FBS_TX_A6(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
FBSRX6 = FBS_RX_A6(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);


EERROR = [max(max(abs(abs(FBSTX0.V) - abs(FBSTX1.V)))) max(max(abs(abs(FBSRX0.V) - abs(FBSRX1.V))));
    max(max(abs(abs(FBSTX0.V) - abs(FBSTX2.V)))) max(max(abs(abs(FBSRX0.V) - abs(FBSRX2.V))));
    max(max(abs(abs(FBSTX0.V) - abs(FBSTX3.V)))) max(max(abs(abs(FBSTX0.V) - abs(FBSRX3.V))));
    max(max(abs(abs(FBSTX0.V) - abs(FBSTX4.V)))) max(max(abs(abs(FBSTX0.V) - abs(FBSRX4.V))));
    max(max(abs(abs(FBSTX0.V) - abs(FBSTX5.V)))) max(max(abs(abs(FBSTX0.V) - abs(FBSRX5.V))));
    max(max(abs(abs(FBSTX0.V) - abs(FBSTX6.V)))) max(max(abs(abs(FBSTX0.V) - abs(FBSRX6.V))))];

DERROR = [max(max(abs(angle(FBSTX0.V) - angle(FBSTX1.V)))) max(max(abs(angle(FBSTX0.V) - angle(FBSRX1.V))));
    max(max(abs(angle(FBSTX0.V) - angle(FBSTX2.V)))) max(max(abs(angle(FBSTX0.V) - angle(FBSRX2.V))));
    max(max(abs(angle(FBSTX0.V) - angle(FBSTX3.V)))) max(max(abs(angle(FBSTX0.V) - angle(FBSRX3.V))));
    max(max(abs(angle(FBSTX0.V) - angle(FBSTX4.V)))) max(max(abs(angle(FBSTX0.V) - angle(FBSRX4.V))));
    max(max(abs(angle(FBSTX0.V) - angle(FBSTX5.V)))) max(max(abs(angle(FBSTX0.V) - angle(FBSRX5.V))));
    max(max(abs(angle(FBSTX0.V) - angle(FBSTX6.V)))) max(max(abs(angle(FBSTX0.V) - angle(FBSRX6.V))))];

VERROR = [max(max(abs(FBSTX0.V - FBSTX1.V)));
    max(max(abs(FBSTX0.V - FBSTX2.V)));
    max(max(abs(FBSTX0.V - FBSTX3.V)));
    max(max(abs(FBSTX0.V - FBSTX4.V)));
    max(max(abs(FBSTX0.V - FBSTX5.V)));
    max(max(abs(FBSTX0.V - FBSTX6.V)))];

PERROR = [max(max(abs(real(FBSTX0.Srx - FBSTX1.Srx))));
    max(max(abs(real(FBSTX0.Srx - FBSTX2.Srx))));
    max(max(abs(real(FBSTX0.Srx - FBSTX3.Srx))));
    max(max(abs(real(FBSTX0.Srx - FBSTX4.Srx))));
    max(max(abs(real(FBSTX0.Srx - FBSTX5.Srx))));
    max(max(abs(real(FBSTX0.Srx - FBSTX6.Srx))))];

QERROR = [max(max(abs(imag(FBSTX0.Srx - FBSTX1.Srx))));
    max(max(abs(imag(FBSTX0.Srx - FBSTX2.Srx))));
    max(max(abs(imag(FBSTX0.Srx - FBSTX3.Srx))));
    max(max(abs(imag(FBSTX0.Srx - FBSTX4.Srx))));
    max(max(abs(imag(FBSTX0.Srx - FBSTX5.Srx))));
    max(max(abs(imag(FBSTX0.Srx - FBSTX6.Srx))))];

SERROR = [max(max(abs(FBSTX0.Srx - FBSTX1.Srx)));
    max(max(abs(FBSTX0.Srx - FBSTX2.Srx)));
    max(max(abs(FBSTX0.Srx - FBSTX3.Srx)));
    max(max(abs(FBSTX0.Srx - FBSTX4.Srx)));
    max(max(abs(FBSTX0.Srx - FBSTX5.Srx)));
    max(max(abs(FBSTX0.Srx - FBSTX6.Srx)))];

% IERROR = [max(max(abs(FBSTX0.I(:,2:end) - FBSTX1.I(:,2:end))));
%     max(max(abs(FBSTX0.I(:,2:end) - FBSTX2.I(:,2:end))));
%     max(max(abs(FBSTX0.I(:,2:end) - FBSTX3.I(:,2:end))));
%     max(max(abs(FBSTX0.I(:,2:end) - FBSTX4.I(:,2:end))));
%     max(max(abs(FBSTX0.I(:,2:end) - FBSTX5.I(:,2:end))));
%     max(max(abs(FBSTX0.I(:,2:end) - FBSTX6.I(:,2:end))))];

%%

% % Plot voltage magnitude across feeder
% figure, box on, hold on   
% for kpath = 1:size(feeder.paths,1)
%     temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
%     plot(temppath,abs(FBSTXbase.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,abs(FBSTXbase.V(2,temppath)),'gx','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,abs(FBSTXbase.V(3,temppath)),'b.','MarkerSize',25,'LineWidth',2)
%     
%     plot(temppath,abs(FBSTXbase.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
%     plot(temppath,abs(FBSTXbase.V(2,temppath)),'g--','MarkerSize',12.5,'LineWidth',1.5)
%     plot(temppath,abs(FBSTXbase.V(3,temppath)),'b--','MarkerSize',25,'LineWidth',1.5)
% end
% plot(1:n,0.95*ones(n),'k--','MarkerSize',12,'LineWidth',2)
% set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',15,'FontWeight','bold')
% % title('Voltage Profile of Base (No Control) Case','FontSize',15,'FontWeight','bold')
% legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
% xlabel('Node','FontSize',15,'FontWeight','bold')
% ylabel('Voltage Magnitude [pu]','FontSize',15,'FontWeight','bold')
% axis([0.5 n+0.5 0.945 1.005])
% pbaspect([1 0.375 1])
% 
% % Plot voltage angle across feeder
% figure, box on, hold on   
% for kpath = 1:size(feeder.paths,1)
%     temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
%     plot(temppath,180/pi*angle(FBSTXbase.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,180/pi*angle(FBSTXbase.V(2,temppath))+120,'gx','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,180/pi*angle(FBSTXbase.V(3,temppath))-120,'b.','MarkerSize',25,'LineWidth',2)
%     
%     plot(temppath,180/pi*angle(FBSTXbase.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
%     plot(temppath,180/pi*angle(FBSTXbase.V(2,temppath))+120,'g--','MarkerSize',12.5,'LineWidth',1.5)
%     plot(temppath,180/pi*angle(FBSTXbase.V(3,temppath))-120,'b--','MarkerSize',25,'LineWidth',1.5)
% end
% set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',15,'FontWeight','bold')
% % title('Feeder Voltage Angle of Base (No Control) Case','FontSize',15,'FontWeight','bold')
% legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
% xlabel('Node','FontSize',15,'FontWeight','bold')
% ylabel('\theta_{n}^{\phi}','FontSize',15,'FontWeight','bold')
% axis([0.5 n+0.5 -2 2])
% pbaspect([1 0.375 1])

