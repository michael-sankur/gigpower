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

%%

V0 = [1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];

FBStest = FBS_3phase_fun_20160105(feeder,loads,controllers,V0);
FBStest.V(FBStest.V == 0) = NaN;

for ph = 1:3
    
    if ph == 1 abc = 'a'; end
    if ph == 2 abc = 'b'; end
    if ph == 3 abc = 'c'; end
    
    figure, box on, hold on
    for kpath = 1:size(feeder.paths,1)
        temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
        plot(temppath,abs(FBStest.V(ph,temppath)),'r.--','MarkerSize',25,'LineWidth',1)
    end
    plot(1:n,0.95*ones(n),'k--','MarkerSize',15,'LineWidth',2)
    set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',15,'FontWeight','bold')
    title(['Phase ' abc ' Voltage Magnitude'],'FontSize',15,'FontWeight','bold')
    legend({['|V_{' abc ',k}|']},...
        'FontSize',15,'FontWeight','bold','location','southwest')
    xlabel('Node','FontSize',15,'FontWeight','bold')
    ylabel('Voltage Magnitude [pu]','FontSize',15,'FontWeight','bold')
    axis([1 n 0.92 1.02])
    pbaspect([1 0.5 1])
    
end

%% Simulation Parameters

sim.iter = 1;

t = 0:1:60;
sim.t = t;

for k1 = 1:length(t)
    loads.spuspu(:,k1,:) = loads.spu(:,:,1) + 2*0.125*sin(4*t(k1)*pi/t(end))*loads.spu(:,:,1) + 4*0.125*(t(k1)/t(end))*loads.spu(:,:,1);
end

%% Simulation Timesteps

resQP = [];
resBase = [];
resCon = [];

tic

for k1 = 1:length(t)
    
    disp('/////////////////////////')
    disp(['k1 = ' num2str(k1)])
    disp('-------------------------')
    
    loads.spu(:,:,1) = loads.spuspu(:,k1,:);
    
    itercount = 0;
    Lfbs = zeros(3,n); sim.Lfbs = Lfbs;
    Vmagiter = ones(3,n);
    Vmagfbs = zeros(3,n);
    
    while max(max(abs(Vmagfbs - Vmagiter))) > 1e-4
    
        itercount = itercount + 1;
        disp(['Iteration: ' num2str(itercount)])

        Vmagiter = Vmagfbs;

        [X,Y,nvar,cvx_optval,cvx_status] = V_Balance_Solver_Yapprox_Losses_20160105(feeder,loads,controllers,sim);
        disp(['CVX Status: ' cvx_status])

        Xa = X(1:nvar);
        Xb = X(nvar+1:2*nvar);
        Xc = X(2*nvar+1:3*nvar);

        XX = [Xa Xb Xc]';

        Yopt = XX(:,1:n);
        Vmagopt = sqrt(Yopt);

        Popt = XX(:,1*n+1:2*n);
        Qopt = XX(:,2*n+1:3*n);
        Sopt = Popt + j*Qopt;

        uopt = XX(:,3*n+1:4*n);
        vopt = XX(:,4*n+1:5*n);
        wopt = uopt + j*vopt;

        controllers.wpu = wopt;
        FBSiter = FBS_3phase_fun_20160105(feeder,loads,controllers,V0);

        for k2 = 2:n
        %     Lfbs(:,k1) = diag(feeder.FZpu(:,:,k1)*FBSiter.I(:,k1)*FBSiter.I(:,k1)');
            Lfbs(:,k2) = diag(feeder.FZpu(:,:,k2)).*diag(FBSiter.I(:,k2)*FBSiter.I(:,k2)');
        end

        sim.Lfbs = Lfbs;
        Vmagfbs = abs(FBSiter.V);

        abs(Vmagfbs - Vmagiter);
        disp(['Magnitude Error: ' num2str(max(max(abs(Vmagfbs - Vmagiter))))])
        disp('-------------------------')

    end
    
    controllers.ww(:,k1,:) = wopt;
        
    controllers.wpu = zeros(3,n);
    FBSbase = FBS_3phase_fun_20160105(feeder,loads,controllers,V0);
    resBase.V(:,k1,:) = FBSbase.V(:,:);
    resBase.I(:,k1,:) = FBSbase.I;
    resBase.Inode(:,k1,:) = FBSbase.Inode;
    resBase.S(:,k1,:) = FBSbase.S;
    resBase.sV(:,k1,:) = FBSbase.sV;
    
    resBase.S0(k1) = sum(abs(FBSbase.S(:,1)));
    resBase.P0(k1) = sum(real(FBSbase.S(:,1)));
    resBase.Q0(k1) = sum(imag(FBSbase.S(:,1)));

    controllers.wpu = wopt;
    FBScon = FBS_3phase_fun_20160105(feeder,loads,controllers,V0);
    resCon.V(:,k1,:) = FBScon.V;
    resCon.I(:,k1,:) = FBScon.I;
    resCon.Inode(:,k1,:) = FBScon.Inode;
    resCon.S(:,k1,:) = FBScon.S;
    resCon.sV(:,k1,:) = FBScon.sV;
    
    resCon.w(:,k1,:) = wopt;

    resCon.S0(k1) = sum(abs(FBScon.S(:,1)));
    resCon.P0(k1) = sum(real(FBScon.S(:,1)));
    resCon.Q0(k1) = sum(imag(FBScon.S(:,1)));
    
    clear X XX Yqp Vmagqp Pqp Qqp Sqp uqp vqp wqp;
    
end

resBase.Vplot = resBase.V;
resBase.Vplot(resBase.Vplot == 0) = NaN;

resCon.Vplot = resCon.V;
resCon.Vplot(resCon.Vplot == 0) = NaN;

toc

%% Compute and plot feeder voltage magnitude imbalance

Vdifbase = zeros(3,length(t),n);
Vdifcon = zeros(3,length(t),n);
Vdifbase_norm = zeros(length(t),n);
Vdifcon_norm = zeros(length(t),n);
for k1 = 1:length(t)
    for knode = 1:n
        if feeder.PH(1,knode) == 1 && feeder.PH(2,knode) == 1 && feeder.PH(3,knode) == 0
            Vdifbase(1,k1,knode) = abs(resBase.V(1,k1,knode)) - abs(resBase.V(2,k1,knode));
            Vdifcon(1,k1,knode) = abs(resCon.V(1,k1,knode)) - abs(resCon.V(2,k1,knode));
        elseif feeder.PH(1,knode) == 0 && feeder.PH(2,knode) == 1 && feeder.PH(3,knode) == 1
            Vdifbase(2,k1,knode) = abs(resBase.V(2,k1,knode)) - abs(resBase.V(3,k1,knode));
            Vdifcon(2,k1,knode) = abs(resCon.V(2,k1,knode)) - abs(resCon.V(3,k1,knode));
        elseif feeder.PH(1,knode) == 1 && feeder.PH(2,knode) == 0 && feeder.PH(3,knode) == 1
            Vdifbase(3,k1,knode) = abs(resBase.V(1,k1,knode)) - abs(resBase.V(3,k1,knode));
            Vdifcon(3,k1,knode) = abs(resCon.V(1,k1,knode)) - abs(resCon.V(3,k1,knode));
        elseif feeder.PH(1,knode) == 1 && feeder.PH(2,knode) == 1 && feeder.PH(3,knode) == 1
            Vdifbase(1,k1,knode) = abs(resBase.V(1,k1,knode)) - abs(resBase.V(2,k1,knode));
            Vdifbase(2,k1,knode) = abs(resBase.V(2,k1,knode)) - abs(resBase.V(3,k1,knode));
            Vdifbase(3,k1,knode) = abs(resBase.V(1,k1,knode)) - abs(resBase.V(3,k1,knode));
            Vdifcon(1,k1,knode) = abs(resCon.V(1,k1,knode)) - abs(resCon.V(2,k1,knode));
            Vdifcon(2,k1,knode) = abs(resCon.V(2,k1,knode)) - abs(resCon.V(3,k1,knode));
            Vdifcon(3,k1,knode) = abs(resCon.V(1,k1,knode)) - abs(resCon.V(3,k1,knode));
        end
        Vdifbase_2norm(knode,k1) = norm(Vdifbase(:,k1,knode),2);
        Vdifcon_2norm(knode,k1) = norm(Vdifcon(:,k1,knode),2);
    end
end

Jbase = sum(Vdifbase_2norm);
Jcon = sum(Vdifcon_2norm);

figure
plot(t,Jbase,'r.',t,Jcon,'g.','MarkerSize',15)
set(gca,'FontSize',15,'FontWeight','bold')
title(strcat('Voltage Imbalance'),'FontSize',15,'FontWeight','bold')
legend({'Base','Control'},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Time [step]','FontSize',15,'FontWeight','bold')
ylabel('J','FontSize',15,'FontWeight','bold')

%% Plot voltage magnitude for a certain node

kplot = 10;

figure, box on
plot(t,abs(resBase.Vplot(1,:,kplot)),'r*',t,abs(resBase.Vplot(2,:,kplot)),'g*',t,abs(resBase.Vplot(3,:,kplot)),'b*',...
    t,abs(resCon.Vplot(1,:,kplot)),'r.',t,abs(resCon.Vplot(2,:,kplot)),'g.',t,abs(resCon.Vplot(3,:,kplot)),'b.',...
    'MarkerSize',15)
set(gca,'YTick',0.9:0.01:1,'FontSize',15,'FontWeight','bold')
title(strcat('Node',' ',feeder.nodelist(kplot),' Voltage Magnitude'),'FontSize',15,'FontWeight','bold')
% legend({['Base'],['Control']},'FontSize',15,'FontWeight','bold','location','southwest')
xlabel('Time [step]','FontSize',15,'FontWeight','bold')
ylabel(strcat('|V_{\phi,',feeder.nodelist(kplot),'}| [pu]'),'FontSize',15,'FontWeight','bold')
axis([t(1) t(end) 0.93 1.01])

%% Plot DER control

close all

for k1 = 1:n
    if sum(controllers.cn(:,k1)) ~= 0
        
        figure(2*k1-1), box on, hold on
        plot(real(controllers.ww(1,:,k1)),imag(controllers.ww(1,:,k1)),'r.',...
            real(controllers.ww(2,:,k1)),imag(controllers.ww(2,:,k1)),'g.',...
            real(controllers.ww(3,:,k1)),imag(controllers.ww(3,:,k1)),'b.','MarkerSize',15)
        plot(controllers.wmaxpu(1,k1)*cos(0:pi/1000:2*pi),controllers.wmaxpu(1,k1)*sin(0:pi/1000:2*pi),'r--',...
            controllers.wmaxpu(2,k1)*cos(0:pi/1000:2*pi),controllers.wmaxpu(2,k1)*sin(0:pi/1000:2*pi),'g--',...
            controllers.wmaxpu(3,k1)*cos(0:pi/1000:2*pi),controllers.wmaxpu(3,k1)*sin(0:pi/1000:2*pi),'b--',...
            'LineWidth',2)
        set(gca,'FontSize',15,'FontWeight','bold')
        legend({['w_{a,' num2str(k1) '}'],['w_{b,' num2str(k1) '}'],['w_{c,' num2str(k1) '}']},'FontSize',15,'FontWeight','bold')
        xlabel(['P_{' num2str(k1) '}'],'FontSize',15,'FontWeight','bold')
        ylabel(['Q_{' num2str(k1) '}'],'FontSize',15,'FontWeight','bold')
        title({['Inverters at node ', feeder.nodelist{k1}]},'FontSize',15,'FontWeight','bold')
        axis([-max(controllers.wmaxpu(:,k1)) max(controllers.wmaxpu(:,k1)) -max(controllers.wmaxpu(:,k1)) max(controllers.wmaxpu(:,k1))]), axis equal
        hold off
        
        print(['-f' num2str(2*k1-1)],'-depsc',['C:\Users\Michael\Desktop\Reconfiguration\Control\figures\eps\wcirc' num2str(num2str(k1)) '.eps'])
        print(['-f' num2str(2*k1-1)],'-dpng',['C:\Users\Michael\Desktop\Reconfiguration\Control\figures\png\wcirc' num2str(num2str(k1)) '.png'])
                
    end
end