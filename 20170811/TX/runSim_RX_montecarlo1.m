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

loads.spu = 1.25*loads.spu;

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

%% IC test

% FBSOPT0 = FBS_RX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom,0);

% FBSOPT1 = FBS_RX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom,1);


%%

FBS0 = FBS_RX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBS1 = FBS_RX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBS2 = FBS_RX_A2(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBS3 = FBS_RX_A3(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBS4 = FBS_RX_A4(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBS5 = FBS_RX_A5(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

FBS6 = FBS_RX_A6(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);

EERROR = [max(max(abs(abs(FBS0.V) - abs(FBS1.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS2.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS3.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS4.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS5.V))));
    max(max(abs(abs(FBS0.V) - abs(FBS6.V))))];

DERROR = [max(max(abs(angle(FBS0.V) - angle(FBS1.V))));
    max(max(abs(angle(FBS0.V) - angle(FBS2.V))));
    max(max(abs(angle(FBS0.V) - angle(FBS3.V))));
    max(max(abs(angle(FBS0.V) - angle(FBS4.V))));
    max(max(abs(angle(FBS0.V) - angle(FBS5.V))));
    max(max(abs(angle(FBS0.V) - angle(FBS6.V))))];

VERROR = [max(max(abs(FBS0.V - FBS1.V)));
    max(max(abs(FBS0.V - FBS2.V)));
    max(max(abs(FBS0.V - FBS3.V)));
    max(max(abs(FBS0.V - FBS4.V)));
    max(max(abs(FBS0.V - FBS5.V)));
    max(max(abs(FBS0.V - FBS6.V)))];

PERROR = [max(max(abs(real(FBS0.Srx - FBS1.Srx))));
    max(max(abs(real(FBS0.Srx - FBS2.Srx))));
    max(max(abs(real(FBS0.Srx - FBS3.Srx))));
    max(max(abs(real(FBS0.Srx - FBS4.Srx))));
    max(max(abs(real(FBS0.Srx - FBS5.Srx))));
    max(max(abs(real(FBS0.Srx - FBS6.Srx))))];

QERROR = [max(max(abs(imag(FBS0.Srx - FBS1.Srx))));
    max(max(abs(imag(FBS0.Srx - FBS2.Srx))));
    max(max(abs(imag(FBS0.Srx - FBS3.Srx))));
    max(max(abs(imag(FBS0.Srx - FBS4.Srx))));
    max(max(abs(imag(FBS0.Srx - FBS5.Srx))));
    max(max(abs(imag(FBS0.Srx - FBS6.Srx))))];

SERROR = [max(max(abs(FBS0.Srx - FBS1.Srx)));
    max(max(abs(FBS0.Srx - FBS2.Srx)));
    max(max(abs(FBS0.Srx - FBS3.Srx)));
    max(max(abs(FBS0.Srx - FBS4.Srx)));
    max(max(abs(FBS0.Srx - FBS5.Srx)));
    max(max(abs(FBS0.Srx - FBS6.Srx)))];

% IERROR = [max(max(abs(FBS0.I(:,2:end) - FBS1.I(:,2:end))));
%     max(max(abs(FBS0.I(:,2:end) - FBS2.I(:,2:end))));
%     max(max(abs(FBS0.I(:,2:end) - FBS3.I(:,2:end))));
%     max(max(abs(FBS0.I(:,2:end) - FBS4.I(:,2:end))));
%     max(max(abs(FBS0.I(:,2:end) - FBS5.I(:,2:end))));
%     max(max(abs(FBS0.I(:,2:end) - FBS6.I(:,2:end))))];

%%

% % Plot voltage magnitude across feeder
% figure, box on, hold on   
% for kpath = 1:size(feeder.paths,1)
%     temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
%     plot(temppath,abs(FBSbase.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,abs(FBSbase.V(2,temppath)),'gx','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,abs(FBSbase.V(3,temppath)),'b.','MarkerSize',25,'LineWidth',2)
%     
%     plot(temppath,abs(FBSbase.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
%     plot(temppath,abs(FBSbase.V(2,temppath)),'g--','MarkerSize',12.5,'LineWidth',1.5)
%     plot(temppath,abs(FBSbase.V(3,temppath)),'b--','MarkerSize',25,'LineWidth',1.5)
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
%     plot(temppath,180/pi*angle(FBSbase.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,180/pi*angle(FBSbase.V(2,temppath))+120,'gx','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,180/pi*angle(FBSbase.V(3,temppath))-120,'b.','MarkerSize',25,'LineWidth',2)
%     
%     plot(temppath,180/pi*angle(FBSbase.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
%     plot(temppath,180/pi*angle(FBSbase.V(2,temppath))+120,'g--','MarkerSize',12.5,'LineWidth',1.5)
%     plot(temppath,180/pi*angle(FBSbase.V(3,temppath))-120,'b--','MarkerSize',25,'LineWidth',1.5)
% end
% set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',15,'FontWeight','bold')
% % title('Feeder Voltage Angle of Base (No Control) Case','FontSize',15,'FontWeight','bold')
% legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
% xlabel('Node','FontSize',15,'FontWeight','bold')
% ylabel('\theta_{n}^{\phi}','FontSize',15,'FontWeight','bold')
% axis([0.5 n+0.5 -2 2])
% pbaspect([1 0.375 1])

