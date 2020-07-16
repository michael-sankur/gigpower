% Michael Sankur - msankur@berkeley.edu
% 2016.01.05

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder. This simulation uses an iterative method to account for losses in
% the power balance.

clc, clear all, close all

% cd C:\Users\Michael\Desktop\cvx-w64\cvx
% cvx_setup

%% Load feeder

feeder = 13;

fpconf = 'C:\Users\Michael\Desktop\Reconfiguration\';
fn1 = ['feeder' num2str(feeder) '_Base.csv'];
fn2 = ['feeder' num2str(feeder) '_Topology.csv'];
fn3 = ['feeder' num2str(feeder) '_Impedance.csv'];
fn4 = ['feeder' num2str(feeder) '_Loads.csv'];
fn5 = ['feeder' num2str(feeder) '_Controllers.csv'];

[feeder, loads, controllers] = feeder_mapper_function_20160105(fpconf,fn1,fn2,fn3,fn4,fn5);

%% Feeder paramaters

n = feeder.n;

feeder.FZpu = 1.25*feeder.FZpu;

%% Load and capacitor parameters

Arand = 0.05*rand(3,n);
A0 = 0.90*ones(3,n) + Arand;
A1 = 0.10*ones(3,n) - Arand;

% A0 = ones(3,n);
% A1 = 0*ones(3,n);

loads.A0 = A0;
loads.A1 = A1;

loads.cappu = 0*loads.cappu;

clear A0 A1 Arand

%% Control parameters

controllers.wpu = zeros(3,n);

%% Simulation parameters

sim.iter = 1;

%%

V0 = [1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];

% FBStest = FBS_3phase_fun_20160105(feeder,loads,controllers,V0);
% FBStest.V(FBStest.V == 0) = NaN;
% 
% for ph = 1:3
%     
%     if ph == 1 abc = 'a'; end
%     if ph == 2 abc = 'b'; end
%     if ph == 3 abc = 'c'; end
%     
%     figure, box on, hold on
%     for kpath = 1:size(feeder.paths,1)
%         temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
%         plot(temppath,abs(FBStest.V(ph,temppath)),'r.--','MarkerSize',25,'LineWidth',1)
%     end
%     plot(1:n,0.95*ones(n),'k--','MarkerSize',15,'LineWidth',2)
%     set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',15,'FontWeight','bold')
%     title(['Phase ' abc ' Voltage Magnitude'],'FontSize',15,'FontWeight','bold')
%     legend({['|V_{' abc ',k}|']},...
%         'FontSize',15,'FontWeight','bold','location','southwest')
%     xlabel('Node','FontSize',15,'FontWeight','bold')
%     ylabel('Voltage Magnitude [pu]','FontSize',15,'FontWeight','bold')
%     axis([1 n 0.92 1.02])
%     pbaspect([1 0.5 1])
%     
% end

%% Run optimization

itercount = 0;
Lfbs = zeros(3,n); sim.Lfbs = Lfbs;
Vmagiter = ones(3,n);
Vmagfbs = zeros(3,n);

while max(max(abs(Vmagfbs - Vmagiter))) > 1e-9
    
    itercount = itercount + 1;
    disp(['Iteration: ' num2str(itercount)])
        
    Vmagiter = Vmagfbs;
    
    [X,Y,nvar,cvx_optval,cvx_status] = V_Balance_Solver_Yapprox_Losses_20160105(feeder,loads,controllers,sim);
    disp(['CVX Status: ' cvx_status])

    Xa = X(1:nvar);
    Xb = X(nvar+1:2*nvar);
    Xc = X(2*nvar+1:3*nvar);

    XX = [Xa Xb Xc]';

    Yqp = XX(:,1:n);
    Vmagqp = sqrt(Yqp);

    Pqp = XX(:,1*n+1:2*n);
    Qqp = XX(:,2*n+1:3*n);
    Sqp = Pqp + j*Qqp;

    uqp = XX(:,3*n+1:4*n);
    vqp = XX(:,4*n+1:5*n);
    wqp = uqp + j*vqp;

    controllers.wpu = wqp;
    FBSiter = FBS_3phase_fun_20160105(feeder,loads,controllers,V0);
    
    for k1 = 2:n
    %     Lfbs(:,k1) = diag(feeder.FZpu(:,:,k1)*FBSiter.I(:,k1)*FBSiter.I(:,k1)');
        Lfbs(:,k1) = diag(feeder.FZpu(:,:,k1)).*diag(FBSiter.I(:,k1)*FBSiter.I(:,k1)');
    end

    sim.Lfbs = Lfbs;
    Vmagfbs = abs(FBSiter.V);
    
    abs(Vmagfbs - Vmagiter);
    disp(['Magnitude Error: ' num2str(max(max(abs(Vmagfbs - Vmagiter))))])
    disp('-------------------------')
    
end

itercount

controllers.wpu = zeros(3,n);
FBSbase = FBS_3phase_fun_20160105(feeder,loads,controllers,V0);

controllers.wpu = wqp;
FBScon = FBS_3phase_fun_20160105(feeder,loads,controllers,V0);


S0_base = sum(abs(FBSbase.S(:,1)))
S0_con = sum(abs(FBScon.S(:,1)))

P0_base = sum(real(FBSbase.S(:,1)))
P0_con = sum(real(FBScon.S(:,1)))

Q0_base = sum(imag(FBSbase.S(:,1)))
Q0_con = sum(imag(FBScon.S(:,1)))

%% Plot feeder voltage magnitude for each phase

close all

fpfig = 'C:\Users\Michael\Desktop\Reconfiguration\IEEE PES 2016\20151118\figures\';
phaselist = ['a','b','c'];

FBSbase.V(FBSbase.V == 0) = NaN;
FBScon.V(FBScon.V == 0) = NaN;

for ph = 1:3
        
    figure, box on, hold on
    for kpath = 1:size(feeder.paths,1)
        temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
        plot(temppath,abs(FBSbase.V(ph,temppath)),'r*','MarkerSize',12.5,'LineWidth',2)
        plot(temppath,abs(FBScon.V(ph,temppath)),'g.','MarkerSize',20,'LineWidth',2)
        plot(temppath,abs(FBSbase.V(ph,temppath)),'r--','MarkerSize',12.5,'LineWidth',0.5)
        plot(temppath,abs(FBScon.V(ph,temppath)),'g--','MarkerSize',20,'LineWidth',0.5)
    end
    plot(0:n+1,0.95*ones(n+2),'k--','LineWidth',2)
    set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'YTick',0.94:0.01:1,'FontSize',15,'FontWeight','bold')
%     title(['Phase ' abc ' Voltage Magnitude'],'FontSize',15,'FontWeight','bold')
%     legend({['|V_{' abc ',k}| - Base'],['|V_{' abc ',k}| - Control']},...
%         'FontSize',15,'FontWeight','bold','location','southwest')
    legend({['Base'],['Control']},...
        'FontSize',15,'FontWeight','bold','location','southwest')
    xlabel('Node','FontSize',15,'FontWeight','bold')
    ylabel(['|V_{' phaselist(ph) ',k}| [pu]'],'FontSize',15,'FontWeight','bold')
    axis([0.5 n+0.5 0.94 1.0])
%     axis([0.5 n+0.5 -inf inf])
    pbaspect([1 0.5 1]) 
    
    print(['-f' num2str(ph)],'-dpng',[fpfig 'png\' 'V' phaselist(ph) '.png'])
    print(['-f' num2str(ph)],'-depsc',[fpfig 'eps\' 'V' phaselist(ph) '.eps'])
    
end

%% Plot feeder head power

% close all

figure, box on, hold on
plot([1 2 3],[sum(abs(FBSbase.S(:,1))) sum(real(FBSbase.S(:,1))) sum(imag(FBSbase.S(:,1)))],'r*',...        
    'MarkerSize',12.5,'LineWidth',2)
plot([1 2 3],[sum(abs(FBScon.S(:,1))) sum(real(FBScon.S(:,1))) sum(imag(FBScon.S(:,1)))],'g.',...
    'MarkerSize',20,'LineWidth',2)
set(gca,'XTick',[1 2 3],'XTickLabel',{'S_0' 'P_0' 'Q_0'},'YTick',0:0.1:1,'FontSize',15,'FontWeight','bold')
% title('Feeder Head Power','FontSize',15,'FontWeight','bold')
legend({['Base'],['Control']},...
    'FontSize',15,'FontWeight','bold','location','southwest')
% xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('Power [pu]','FontSize',15,'FontWeight','bold')
% axis([0.5 3.5 0.2 1])
axis([0.5 3.5 -inf inf])
pbaspect([1 0.5 1]) 

print('-f4','-dpng',[fpfig 'png\' 'SPQ' '.png'])
print('-f4','-depsc',[fpfig 'eps\' 'SPQ' '.eps'])

%% Plot feeder voltage magntiude on all phases

% close all

figure, box on, hold on
plot(1:n,abs(FBSbase.V(1,:)),'r*','MarkerSize',12.5,'LineWidth',2)
plot(1:n,abs(FBSbase.V(2,:)),'gd','MarkerSize',12.5,'LineWidth',2)
plot(1:n,abs(FBSbase.V(3,:)),'b.','MarkerSize',20,'LineWidth',2)
plot(0:n+1,0.95*ones(n+2),'k--','LineWidth',2)
set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'YTick',0.94:0.01:1,'FontSize',15,'FontWeight','bold')
% title(['Base Scenario Voltage Magnitudes'],'FontSize',15,'FontWeight','bold')
legend({'|V_{a,k}|','|V_{b,k}|','|V_{c,k}|'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('|V_{\phi,k}| [pu]','FontSize',15,'FontWeight','bold')
axis([0.5 n+0.5 0.94 1.0])
% axis([0.5 n+0.5 -inf inf])
pbaspect([1 0.5 1]) 

figure, box on, hold on
plot(1:n,abs(FBScon.V(1,:)),'r*','MarkerSize',12.5,'LineWidth',2)
plot(1:n,abs(FBScon.V(2,:)),'gd','MarkerSize',12.5,'LineWidth',2)
plot(1:n,abs(FBScon.V(3,:)),'b.','MarkerSize',20,'LineWidth',2)
plot(0:n+1,0.95*ones(n+2),'k--','LineWidth',2)
set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'YTick',0.94:0.01:1,'FontSize',15,'FontWeight','bold')
% title(['Control Scenario Voltage Magnitudes'],'FontSize',15,'FontWeight','bold')
legend({'|V_{a,k}|','|V_{b,k}|','|V_{c,k}|'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('|V_{\phi,k}| [pu]','FontSize',15,'FontWeight','bold')
axis([0.5 n+0.5 0.94 1.0])
% axis([0.5 n+0.5 -inf inf])
pbaspect([1 0.5 1]) 

print('-f5','-dpng',[fpfig 'png\' 'Vbase' '.png'])
print('-f5','-depsc',[fpfig 'eps\' 'Vbase' '.eps'])

print('-f6','-dpng',[fpfig 'png\' 'Vbase' '.png'])
print('-f6','-depsc',[fpfig 'eps\' 'Vbase' '.eps'])

%% Plot feeder voltage magnitude imbalance

% close all

% rows - a - b ; b - c; a - c
% columns - node

Vdifbase = zeros(3,n);
Vdifcon = zeros(3,n);
Vdifbase_norm = zeros(n);
Vdifcon_norm = zeros(n);
for knode = 1:n
    if feeder.PH(1,knode) == 1 && feeder.PH(2,knode) == 1 && feeder.PH(3,knode) == 0
        Vdifbase(1,knode) = abs(FBSbase.V(1,knode)) - abs(FBSbase.V(2,knode));
        Vdifcon(1,knode) = abs(FBScon.V(1,knode)) - abs(FBScon.V(2,knode));
    elseif feeder.PH(1,knode) == 0 && feeder.PH(2,knode) == 1 && feeder.PH(3,knode) == 1
        Vdifbase(2,knode) = abs(FBSbase.V(2,knode)) - abs(FBSbase.V(3,knode));
        Vdifcon(2,knode) = abs(FBScon.V(2,knode)) - abs(FBScon.V(3,knode));
    elseif feeder.PH(1,knode) == 1 && feeder.PH(2,knode) == 0 && feeder.PH(3,knode) == 1
        Vdifbase(3,knode) = abs(FBSbase.V(1,knode)) - abs(FBSbase.V(3,knode));
        Vdifcon(3,knode) = abs(FBScon.V(1,knode)) - abs(FBScon.V(3,knode));
    elseif feeder.PH(1,knode) == 1 && feeder.PH(2,knode) == 1 && feeder.PH(3,knode) == 1
        Vdifbase(1,knode) = abs(FBSbase.V(1,knode)) - abs(FBSbase.V(2,knode));
        Vdifbase(2,knode) = abs(FBSbase.V(2,knode)) - abs(FBSbase.V(3,knode));
        Vdifbase(3,knode) = abs(FBSbase.V(1,knode)) - abs(FBSbase.V(3,knode));
        Vdifcon(1,knode) = abs(FBScon.V(1,knode)) - abs(FBScon.V(2,knode));
        Vdifcon(2,knode) = abs(FBScon.V(2,knode)) - abs(FBScon.V(3,knode));
        Vdifcon(3,knode) = abs(FBScon.V(1,knode)) - abs(FBScon.V(3,knode));
    end
    Vdifbase_2norm(knode) = norm(Vdifbase(:,knode),2);
    Vdifcon_2norm(knode) = norm(Vdifcon(:,knode),2);
end

% Vdifbase_2norm(Vdifbase_2norm == 0) = NaN;
% Vdifcon_2norm(Vdifcon_2norm == 0) = NaN;

Jbase = sum(Vdifbase_2norm)
Jcon = sum(Vdifcon_2norm)

figure, box on
semilogy(1:n,Vdifbase_2norm,'r*','MarkerSize',12.5,'LineWidth',2), hold on
semilogy(1:n,Vdifcon_2norm,'g.','MarkerSize',20,'LineWidth',2), hold on
set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'YTick',logspace(-4,-1,4),'FontSize',15,'FontWeight','bold')
% title(['Voltage Magnitude Imbalance'],'FontSize',15,'FontWeight','bold')
legend({'Base','Control'},'FontSize',15,'FontWeight','bold','location','east')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('J_{k}','FontSize',15,'FontWeight','bold')
axis([0.5 n+0.5 1e-4 1e-1])
pbaspect([1 0.5 1]) 

print('-f7','-dpng',[ fpfig 'png\' 'Vbalance' '.png'])
print('-f7','-depsc',[ fpfig 'eps\' 'Vbalance' '.eps'])

%% Plot control effort

% close all

figure, box on, hold on
plot(1:length(controllers.cnodes),imag(wqp(1,controllers.cnodes)),'r*','MarkerSize',12.5,'LineWidth',2)
plot(1:length(controllers.cnodes),imag(wqp(2,controllers.cnodes)),'gd','MarkerSize',12.5,'LineWidth',2)
plot(1:length(controllers.cnodes),imag(wqp(3,controllers.cnodes)),'b.','MarkerSize',20,'LineWidth',2)
set(gca,'XTick',1:length(controllers.cnodes),'XTickLabel',feeder.nodelist(controllers.cnodes),'YTick',-0.1:0.025:0.1,'FontSize',15,'FontWeight','bold')
% title(['DER Control Input'],'FontSize',15,'FontWeight','bold')
legend({'u_{a,k}','u_{b,k}','u_{c,k}'},'FontSize',15,'FontWeight','bold','location','northwest')
xlabel('Node','FontSize',15,'FontWeight','bold')
ylabel('[VAR pu]','FontSize',15,'FontWeight','bold')
axis([0 length(controllers.cnodes)+1 -max(max(controllers.wmaxpu)) max(max(controllers.wmaxpu))])
pbaspect([1 0.5 1])

print('-f8','-dpng',[ fpfig 'png\' 'unode' '.png'])
print('-f8','-depsc',[ fpfig 'eps\' 'unode' '.eps'])
