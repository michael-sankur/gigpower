% Michael Sankur - msankur@lbl.gov
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

feeder.name = 'ieee_13node';

% Read from csv file, input first 4 columns as strings (general case)
% Change file path as needed
fpconf = ['C:\Users\Michael\Desktop\Reconfiguration\Feeders\' feeder.name '\'];
fn1 = [feeder.name '_Base.csv.'];
fn2 = [feeder.name '_Topology.csv.'];
fn3 = [feeder.name '_Impedance.csv.'];
fn4 = [feeder.name '_Loads.csv.'];
fn5 = [feeder.name '_Controllers.csv.'];

[feeder1, loads1, controllers1] = feeder_mapper_ieee_function_20160603(fpconf,fn1,fn2,fn3,fn4,fn5);

%% Feeder paramaters

n1 = feeder1.n;

% feeder impedance matrix of dimension 3x3xn [pu]
% feeder.FZpu = 2.0*feeder.FZpu;
feeder1.FZpu = 1.0*feeder1.FZpu;

%% Load and capacitor parameters

% Parameters for voltage dependent loads
% s(y_m,k) = A0(m,k) + A1(m,n)*y_m,k
% Arand = 0.05*rand(3,n);
% A0 = 0.90*ones(3,n) + Arand;
% A1 = 0.10*ones(3,n) - Arand;
% 
% A0 = ones(3,n);
% A1 = 0*ones(3,n);
% 
% loads.A0 = A0;
% loads.A1 = A1;

% Modify ZIP model paramaters
for k1 = 1:n1
    if strmatch(loads1.type{k1},'CI','exact')
        loads1.APQ(:,k1) = 0.9*ones(3,1);
        loads1.AZ(:,k1) = 0.1*ones(3,1);
    elseif strmatch(loads1.type{k1},'CPQ','exact')
        loads1.APQ(:,k1) = ones(3,1);
        loads1.AZ(:,k1) = zeros(3,1);
    elseif strmatch(loads1.type{k1},'CZ','exact')
        loads1.APQ(:,k1) = zeros(3,1);
        loads1.AZ(:,k1) = ones(3,1);
    end
end

% Modify ZIP model paramaters again
loads1.APQ = 1*ones(3,n1).*feeder1.PH;
loads1.AI = 0*ones(3,n1).*feeder1.PH;
loads1.AZ = 0*ones(3,n1).*feeder1.PH;

% feeder loads s (wye connected) at nodes [pu]
loads1.s = 1.0*loads1.s;
loads1.spu = 1.0*loads1.spu;

% feeder capacitance at nodes [pu]
loads1.cappu = 1*loads1.cappu;

clear A0 A1 Arand

%% Control parameters

% current control action of DER [pu]
controllers1.wpu = zeros(3,n1);

% previous node DER control [pu]
controllers1.wk1 = zeros(3,n1);

% max control bounds on DER control [pu]
controllers1.wmax = 1.5*controllers1.wmax;
controllers1.wmaxpu = 1.5*controllers1.wmaxpu;

%% Simulation parameters

V0 = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim1.V0 = V0;

sim1.Vfbs = [];
sim1.Lfbs = [];
sim1.Hfbs = [];

%% FBS Old

FBS1 = FBS_3phase_fun_20160603(feeder1,loads1,controllers1,V0)

% FBS1.V(FBSold.V == 0) = NaN;

%% Load feeder

% name = '6node_fullphase_test';
name = 'ieee_13node_new';

fp = ['C:\Users\Michael\Desktop\reconfiguration\Feeders\'];
fn = [name '.txt'];

[feeder2, nodes2, lines2, configs2, loads2, caps2, controllers2] = feeder_mapper_TX_function_20170215(name, fp, fn);


%% Feeder paramaters

nnode = nodes2.nnode;
nline = lines2.nline;

%% Load parameters

loads2.aPQ = ones(3,nnode).*nodes2.PH;
loads2.aI = zeros(3,nnode);
loads2.aZ = zeros(3,nnode);

loads2.spu = loads2.spu;

%% Controller parameters

controllers2.wmaxpu = 1*controllers2.wmaxpu;

controllers2.wpu = zeros(3,nnode);

%% FBS New

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim2.Vnom = Vnom;

sim2.Vfbs = [];
sim2.Lfbs = [];
sim2.Hfbs = [];

sim.rho = 0.1;

FBS2 = FBS_RX(feeder2,nodes2,lines2,configs2,loads2,caps2,controllers2,Vnom)

% for k1 = 2:lines2.nline
%     lines2.FZpu(:,:,k1) = feeder1.FZpu(:,:,1:k1-1);
% end
lines2.FZpu(:,:,2:end) = feeder1.FZpu(:,:,2:end);


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

