% Michael Sankur - msankur@berkeley.edu
% 2016.04.06

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

% cd C:\Users\Michael\Desktop\cvx-w64\cvx
% cvx_setup

% addpath('C:\Users\Michael\Desktop\columnlegend\')

%% Load feeder

feeder.name = 'feeder13_ieee_cdc';

% Read from csv file, input first 4 columns as strings (general case)
% Change file path as needed
fpconf = ['C:\Users\Michael\Desktop\Reconfiguration\' feeder.name '\'];
fn1 = [feeder.name '_Base.csv.'];
fn2 = [feeder.name '_Topology.csv.'];
fn3 = [feeder.name '_Impedance.csv.'];
fn4 = [feeder.name '_Loads.csv.'];
fn5 = [feeder.name '_Controllers.csv.'];

[feeder, loads, controllers] = feeder_mapper_function_20160603(fpconf,fn1,fn2,fn3,fn4,fn5);

%% Feeder paramaters

n = feeder.n;

% feeder impedance matrix of dimension 3x3xn [pu]
% feeder.FZpu = 2.0*feeder.FZpu;
feeder.FZpu = 1.0*feeder.FZpu;

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

for k1 = 1:n    
    if strmatch(loads.type{k1},'CI','exact')
        loads.APQ(:,k1) = 0.9*ones(3,1);
        loads.AZ(:,k1) = 0.1*ones(3,1);
    elseif strmatch(loads.type{k1},'CPQ','exact')
        loads.APQ(:,k1) = ones(3,1);
        loads.AZ(:,k1) = zeros(3,1);
    elseif strmatch(loads.type{k1},'CZ','exact')
        loads.APQ(:,k1) = zeros(3,1);
        loads.AZ(:,k1) = ones(3,1);
    end
end

loads.APQ = 0.9*ones(3,n);
loads.AI = 0*ones(3,n);
loads.AZ = 0.1*ones(3,n);

% feeder loads s (wye connected) at nodes [pu]
% loads.spu = 4*1.25*loads.spu;
loads.s = 1.25*loads.s;
loads.spu = 1.25*loads.spu;

% loads.spu(3,2) = 0.25*loads.spu(3,2);

% feeder capacitance at nodes [pu]
loads.cappu = 0*loads.cappu;

clear A0 A1 Arand

%% Control parameters

% control of DER at nodes [pu]
controllers.wpu = zeros(3,n);

% previous node DER control [pu]
controllers.wk1 = zeros(3,n);

% max control bounds on DER control [pu]
% can be used for square and magnitude bounds
controllers.wmax = 1*controllers.wmax;
controllers.wmaxpu = 1*controllers.wmaxpu;

%% Simulation parameters

V0 = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim.V0 = V0;

sim.Vfbs = [];
sim.Lfbs = [];
sim.Hfbs = [];

tn = 7;

mag_ref = NaN*ones(3,n);
mag_ref(:,tn) = [1 1 1]';
sim.mag_ref = mag_ref;

delta_ref = NaN*ones(3,n);
delta_ref(:,tn) = [0 -120 120]';
sim.delta_ref = delta_ref;

delta0 = [0 -120 120];
% sim.delta0 = delta0;

Vref = mag_ref.*exp(j*pi/180*delta_ref);

%% Plot feeder state without control

FBSbase = FBS_3phase_fun_20160603(feeder,loads,controllers,V0);

%% Obtain Dist3Flow solution

disp('----- Dist3Flow Formulation -----')

sim.iter = 0;
sim.ramp = 0;
optsol = Solver_LinDist3Flow_Tracking_20160603(feeder,loads,controllers,sim);

disp(['CVX Status: ' optsol.cvx_status])
disp(['Dist3Flow y Cost: ' num2str(optsol.Zy)])
disp(['Dist3Flow theta Cost: ' num2str(optsol.Ztheta)])
disp(['Dist3Flow w Cost: ' num2str(optsol.Zw)])

nvar = optsol.nvar;
Xa = optsol.X(1:nvar);
Xb = optsol.X(nvar+1:2*nvar);
Xc = optsol.X(2*nvar+1:3*nvar);

XX = [Xa Xb Xc]';

Yqp = XX(:,1:n);
Vmagapp = XX(:,1*n+1:2*n);
Vmagqp = sqrt(Yqp);

Dqp = XX(:,2*n+1:3*n);

Pqp = XX(:,3*n+1:4*n);
Qqp = XX(:,4*n+1:5*n);
Sqp = Pqp + j*Qqp;

uqp = XX(:,5*n+1:6*n);
vqp = XX(:,6*n+1:7*n);
wqp = uqp + j*vqp;

controllers.wpu = wqp;
FBSd3f = FBS_3phase_fun_20160603(feeder,loads,controllers,V0);

disp('/////////////////////////')
disp(' ')

%% Plot voltage magnitude for base case

FBSbase.V(FBSbase.V == 0) = NaN;

figure, box on, hold on   
for kpath = 1:size(feeder.paths,1)
    temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
    plot(temppath,abs(FBSbase.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
    plot(temppath,abs(FBSbase.V(2,temppath)),'gx','MarkerSize',12.5,'LineWidth',2)
    plot(temppath,abs(FBSbase.V(3,temppath)),'b.','MarkerSize',25,'LineWidth',2)
    
    plot(temppath,abs(FBSbase.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
    plot(temppath,abs(FBSbase.V(2,temppath)),'g--','MarkerSize',12.5,'LineWidth',1.5)
    plot(temppath,abs(FBSbase.V(3,temppath)),'b--','MarkerSize',25,'LineWidth',1.5)
end
plot(1:n,0.95*ones(n),'k--','MarkerSize',12,'LineWidth',2)
set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',15,'FontWeight','bold')
% title('Voltage Profile of Base (No Control) Case','FontSize',15,'FontWeight','bold')
legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('Voltage Magnitude [pu]','FontSize',15,'FontWeight','bold')
axis([0.5 n+0.5 0.945 1.005])
pbaspect([1 0.375 1])

% print('-f1','-dpng','C:\Users\Michael\Desktop\temp\37node_base.png')
% print('-f1','-depsc','C:\Users\Michael\Desktop\temp\37node_base.eps')

%% Plot voltage magnitude from FBS with SDP control

FBSd3f.V(feeder.PH == 0) = NaN;

figure, box on, hold on   
for kpath = 1:size(feeder.paths,1)
    temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
    plot(temppath,abs(FBSd3f.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
    plot(temppath,abs(FBSd3f.V(2,temppath)),'gx','MarkerSize',12.5,'LineWidth',2)
    plot(temppath,abs(FBSd3f.V(3,temppath)),'b.','MarkerSize',25,'LineWidth',2)
    
    plot(temppath,abs(FBSd3f.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
    plot(temppath,abs(FBSd3f.V(2,temppath)),'g--','MarkerSize',12.5,'LineWidth',1.5)
    plot(temppath,abs(FBSd3f.V(3,temppath)),'b--','MarkerSize',25,'LineWidth',1.5)
end
for k1 = 2:n
    plot(k1,mag_ref(1,k1),'rs','MarkerSize',15,'LineWidth',1.5)
    plot(k1,mag_ref(2,k1),'gs','MarkerSize',15,'LineWidth',1.5)
    plot(k1,mag_ref(3,k1),'bs','MarkerSize',15,'LineWidth',1.5)
end
plot(1:n,0.95*ones(n),'k--','MarkerSize',12,'LineWidth',2)
set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',15,'FontWeight','bold')
% title('Feeder voltage magnitude using D3F control','FontSize',15,'FontWeight','bold')
legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('| V_{m}^{\phi} | [pu]','FontSize',15,'FontWeight','bold')
axis([0.5 n+0.5 0.95 1.005])
pbaspect([1 0.375 1])

% print('-f3','-dpng','C:\Users\Michael\Desktop\temp\sdpcomp_mag_d3f.png')
% print('-f3','-depsc','C:\Users\Michael\Desktop\temp\sdpcomp_mag_d3f.eps')

%%

FBSbase.V(FBSbase.V == 0) = NaN;

figure, box on, hold on   
for kpath = 1:size(feeder.paths,1)
    temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
    plot(temppath,180/pi*angle(FBSbase.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
    plot(temppath,180/pi*angle(FBSbase.V(2,temppath))+120,'gx','MarkerSize',12.5,'LineWidth',2)
    plot(temppath,180/pi*angle(FBSbase.V(3,temppath))-120,'b.','MarkerSize',25,'LineWidth',2)
    
    plot(temppath,180/pi*angle(FBSbase.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
    plot(temppath,180/pi*angle(FBSbase.V(2,temppath))+120,'g--','MarkerSize',12.5,'LineWidth',1.5)
    plot(temppath,180/pi*angle(FBSbase.V(3,temppath))-120,'b--','MarkerSize',25,'LineWidth',1.5)
end
set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',15,'FontWeight','bold')
% title('Feeder Voltage Angle of Base (No Control) Case','FontSize',15,'FontWeight','bold')
legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('\theta_{m}^{\phi}','FontSize',15,'FontWeight','bold')
axis([0.5 n+0.5 -2 2])
pbaspect([1 0.375 1])

% print('-f4','-dpng','C:\Users\Michael\Desktop\temp\sdpcomp_ang_base.png')
% print('-f4','-depsc','C:\Users\Michael\Desktop\temp\sdpcomp_ang_base.eps')

%%

FBSd3f.V(feeder.PH == 0) = NaN;

figure, box on, hold on
for kpath = 1:size(feeder.paths,1)
    temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
    plot(temppath,180/pi*angle(FBSd3f.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
    plot(temppath,180/pi*angle(FBSd3f.V(2,temppath))+120,'gx','MarkerSize',12.5,'LineWidth',2)
    plot(temppath,180/pi*angle(FBSd3f.V(3,temppath))-120,'b.','MarkerSize',25,'LineWidth',2)
    
    plot(temppath,180/pi*angle(FBSd3f.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
    plot(temppath,180/pi*angle(FBSd3f.V(2,temppath))+120,'g--','MarkerSize',12.5,'LineWidth',1.5)
    plot(temppath,180/pi*angle(FBSd3f.V(3,temppath))-120,'b--','MarkerSize',25,'LineWidth',1.5)
end
for k1 = 2:n
    plot(k1,delta_ref(1,k1)-delta0(1),'rd','MarkerSize',15,'LineWidth',1.5)
    plot(k1,delta_ref(2,k1)-delta0(2),'gd','MarkerSize',15,'LineWidth',1.5)
    plot(k1,delta_ref(3,k1)-delta0(3),'bd','MarkerSize',15,'LineWidth',1.5)
end
set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',15,'FontWeight','bold')
% title('Voltage Angle using D3F Control','FontSize',15,'FontWeight','bold')
legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('\theta_{0n}^{\phi}','FontSize',15,'FontWeight','bold')
axis([0.5 n+0.5 -2 2])
pbaspect([1 0.375 1])

% print('-f6','-dpng','C:\Users\Michael\Desktop\temp\sdpcomp_ang_d3f.png')
% print('-f6','-depsc','C:\Users\Michael\Desktop\temp\sdpcomp_ang_d3f.eps')

%%

kpath = tn;
temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
len = length(temppath);

figure, box on, hold on
plot(1:len,abs(FBSbase.V(1,temppath)),'r+','MarkerSize',15,'LineWidth',2)
plot(1:len,abs(FBSbase.V(2,temppath)),'gx','MarkerSize',15,'LineWidth',2)
plot(1:len,abs(FBSbase.V(3,temppath)),'b.','MarkerSize',25,'LineWidth',2)

plot(1:len,abs(FBSbase.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
plot(1:len,abs(FBSbase.V(2,temppath)),'g--','MarkerSize',12.5,'LineWidth',1.5)
plot(1:len,abs(FBSbase.V(3,temppath)),'b--','MarkerSize',25,'LineWidth',1.5)
plot(1:n,0.95*ones(n),'k--','MarkerSize',12,'LineWidth',2)
set(gca,'XTick',1:len,'XTickLabel',feeder.nodelist(temppath),'FontSize',15,'FontWeight','bold')
% title('Voltage Profile of Base (No Control) Case','FontSize',15,'FontWeight','bold')
legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('Voltage Magnitude [pu]','FontSize',15,'FontWeight','bold')
axis([0.5 len+0.5 0.945 1.005])
pbaspect([1 0.375 1])

print('-f5','-dpng',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'path_mag_base.png'])
print('-f5','-depsc',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'path_mag_base.eps'])

%%

figure, box on, hold on
plot(1:len,abs(FBSd3f.V(1,temppath)),'r+','MarkerSize',15,'LineWidth',2)
plot(1:len,abs(FBSd3f.V(2,temppath)),'gx','MarkerSize',15,'LineWidth',2)
plot(1:len,abs(FBSd3f.V(3,temppath)),'b.','MarkerSize',25,'LineWidth',2)

plot(len,mag_ref(1,temppath(end)),'ks','MarkerSize',15,'LineWidth',2)

plot(1:len,abs(FBSd3f.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
plot(1:len,abs(FBSd3f.V(2,temppath)),'g--','MarkerSize',12.5,'LineWidth',1.5)
plot(1:len,abs(FBSd3f.V(3,temppath)),'b--','MarkerSize',25,'LineWidth',1.5)
plot(1:n,0.95*ones(n),'k--','MarkerSize',12,'LineWidth',2)
set(gca,'XTick',1:len,'XTickLabel',feeder.nodelist(temppath),'FontSize',15,'FontWeight','bold')
% title('Voltage Profile of Base (No Control) Case','FontSize',15,'FontWeight','bold')
legend({'a','b','c','Magnitude Reference'},'FontSize',15,'FontWeight','bold','location','southeast')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('Voltage Magnitude [pu]','FontSize',15,'FontWeight','bold')
axis([0.5 len+0.5 0.985 1.005])
pbaspect([1 0.375 1])

print('-f6','-dpng',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'path_mag_control.png'])
print('-f6','-depsc',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'path_mag_control.eps'])

%%

figure, box on, hold on
temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
plot(1:len,180/pi*angle(FBSbase.V(1,temppath)),'r+','MarkerSize',15,'LineWidth',2)
plot(1:len,180/pi*angle(FBSbase.V(2,temppath))+120,'gx','MarkerSize',15,'LineWidth',2)
plot(1:len,180/pi*angle(FBSbase.V(3,temppath))-120,'b.','MarkerSize',25,'LineWidth',2)

plot(1:len,180/pi*angle(FBSbase.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
plot(1:len,180/pi*angle(FBSbase.V(2,temppath))+120,'g--','MarkerSize',12.5,'LineWidth',1.5)
plot(1:len,180/pi*angle(FBSbase.V(3,temppath))-120,'b--','MarkerSize',25,'LineWidth',1.5)
set(gca,'XTick',1:len,'XTickLabel',feeder.nodelist(temppath),'YTick',-1:0.1:1,'FontSize',15,'FontWeight','bold')
% title('Feeder Voltage Angle of Base (No Control) Case','FontSize',15,'FontWeight','bold')
legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('Voltage Angle [\circ]','FontSize',15,'FontWeight','bold')
axis([0.5 len+0.5 -0.9 0.1])
pbaspect([1 0.375 1])

print('-f7','-dpng',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'path_angle_base.png'])
print('-f7','-depsc',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'path_angle_base.eps'])

%%

figure, box on, hold on
temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
plot(1:len,180/pi*angle(FBSd3f.V(1,temppath)),'r+','MarkerSize',15,'LineWidth',2)
plot(1:len,180/pi*angle(FBSd3f.V(2,temppath))+120,'gx','MarkerSize',15,'LineWidth',2)
plot(1:len,180/pi*angle(FBSd3f.V(3,temppath))-120,'b.','MarkerSize',25,'LineWidth',2)

plot(len,delta_ref(1,temppath(end)),'ks','MarkerSize',15,'LineWidth',2)

plot(1:len,180/pi*angle(FBSd3f.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
plot(1:len,180/pi*angle(FBSd3f.V(2,temppath))+120,'g--','MarkerSize',12.5,'LineWidth',1.5)
plot(1:len,180/pi*angle(FBSd3f.V(3,temppath))-120,'b--','MarkerSize',25,'LineWidth',1.5)

set(gca,'XTick',1:len,'XTickLabel',feeder.nodelist(temppath),'YTick',-1:0.05:1,'FontSize',15,'FontWeight','bold')
% title('Feeder Voltage Angle of Base (No Control) Case','FontSize',15,'FontWeight','bold')
legend({'a','b','c','Angle Reference'},'FontSize',15,'FontWeight','bold','location','southeast')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('Voltage Angle [\circ]','FontSize',15,'FontWeight','bold')
% axis([0.5 len+0.5 -0.15 0.05])
pbaspect([1 0.375 1])

print('-f8','-dpng',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'path_angle_control.png'])
print('-f8','-depsc',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'path_angle_control.eps'])

%%

figure, box on, hold on, axis equal
plot(real(FBSbase.V(1,tn)),imag(FBSbase.V(1,tn)),'r*', ...
    real(FBSbase.V(2,tn)*exp(j*120*pi/180)),imag(FBSbase.V(2,tn)*exp(j*120*pi/180)),'g*', ...
    real(FBSbase.V(3,tn)*exp(j*-120*pi/180)),imag(FBSbase.V(3,tn)*exp(j*-120*pi/180)),'b*', ...
    'MarkerSize',15,'LineWidth',2)
plot(real(FBSd3f.V(1,tn)),imag(FBSd3f.V(1,tn)),'r.', ...
    real(FBSd3f.V(2,tn)*exp(j*120*pi/180)),imag(FBSd3f.V(2,tn)*exp(j*120*pi/180)),'g.', ...
    real(FBSd3f.V(3,tn)*exp(j*-120*pi/180)),imag(FBSd3f.V(3,tn)*exp(j*-120*pi/180)),'b.', ...
    'MarkerSize',25,'LineWidth',2)
 plot(real(mag_ref(1,tn)),imag(mag_ref(1,tn)),'kd', ...
    'MarkerSize',12.5,'LineWidth',2)
set(gca,'FontSize',15,'FontWeight','bold')
% title('Feeder Voltage Angle of Base (No Control) Case','FontSize',15,'FontWeight','bold')
% legend({'a - base','b - base','c - base','a - control','b - control','c - control','reference'}, ...
%     'FontSize',15,'FontWeight','bold','location','northwest')
xlabel('Real','FontSize',15,'FontWeight','bold')
ylabel('Imag','FontSize',15,'FontWeight','bold')
axis([0.945 1.005 -0.015 0.015])
% pbaspect([1 0.375 1])

legend_str{1} = 'a - base';
legend_str{2} = 'b - base';
legend_str{3} = 'c - base';
legend_str{4} = 'a - control';
legend_str{5} = 'b - control';
legend_str{6} = 'c - control';
% legend_str{7} = 'reference';
columnlegend(2,legend_str,'boxon','Location','NorthWest','FontSize',15,'FontWeight','bold');

print('-f9','-dpng',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'node_real_imag.png'])
print('-f9','-depsc',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'node_real_imag.eps'])

%%

figure, box on, axis equal
plot(real(mag_ref(1,tn)),imag(mag_ref(1,tn)),'k*',real(FBSd3f.V(1,tn)),imag(FBSd3f.V(1,tn)),'r.', ...
    real(FBSd3f.V(2,tn)*exp(j*120*pi/180)),imag(FBSd3f.V(2,tn)*exp(j*120*pi/180)),'g.', ...
    real(FBSd3f.V(3,tn)*exp(j*-120*pi/180)),imag(FBSd3f.V(3,tn)*exp(j*-120*pi/180)),'b.', ...
    0.995*cosd((-0.1:0.001:0.1)),0.995*sind((-0.1:0.001:0.1)),'k--', ...
    1.005*cosd((-0.1:0.001:0.1)),1.005*sind((-0.1:0.001:0.1)),'k--', ...
    (0.995:0.0001:1.005)*cosd(-0.1),(0.995:0.0001:1.005)*sind(-0.1),'k--', ...
    (0.995:0.0001:1.005)*cosd(0.1),(0.995:0.0001:1.005)*sind(0.1),'k--', ...
    'MarkerSize',15,'LineWidth',1.5)
set(gca,'FontSize',15,'FontWeight','bold')
legend({'Target','V_{3}^{a}','V_{3}^{b}','bounds'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Real','FontSize',15,'FontWeight','bold')
ylabel('Imag','FontSize',15,'FontWeight','bold')
% axis([0.5 n+0.5 -1 0.1])
% pbaspect([1 0.375 1])

% print('-f6','-dpng','C:\Users\Michael\Desktop\temp\sdpcomp_ang_d3f.png')
% print('-f6','-depsc','C:\Users\Michael\Desktop\temp\sdpcomp_ang_d3f.eps')

%%

figure, box on
plot(mag_ref(1,tn) - abs(FBSd3f.V(1,tn)),delta_ref(1,tn) - 180/pi*angle(FBSd3f.V(1,tn)),'r*', ...
    mag_ref(2,tn) - abs(FBSd3f.V(2,tn)),delta_ref(2,tn) - 180/pi*angle(FBSd3f.V(2,tn)),'g*', ...
    mag_ref(3,tn) - abs(FBSd3f.V(3,tn)),delta_ref(3,tn) - 180/pi*angle(FBSd3f.V(3,tn)),'b*', ...
    'MarkerSize',20,'LineWidth',1.5)
set(gca,'XTick',-0.006:0.001:0.006,'YTick',-0.010:0.001:0.010,'FontSize',15,'FontWeight','bold')
legend({'a','b','c'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Magnitude Error','FontSize',15,'FontWeight','bold')
ylabel('Angle Error','FontSize',15,'FontWeight','bold')
axis([0 0.003 -0.003 0.009])
pbaspect([1 0.5 1])

print('-f11','-dpng',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'node_error_mag_angle.png'])
print('-f11','-depsc',['C:\Users\Michael\Desktop\temp\track_' char(feeder.nodelist(tn)) 'node_error_mag_angle.eps'])
