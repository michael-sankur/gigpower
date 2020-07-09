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

name = '6node_fullphase_test';
% name = 'ieee_13node_new';


fp = ['C:\Users\Michael\Desktop\reconfiguration\Feeders\'];
fn = [name '.txt'];

[feeder, nodes, lines, configs, loads, caps, controllers] = feeder_mapper_TX_function_20170215(name, fp, fn);


%% Feeder paramaters

nnode = nodes.nnode;
nline = lines.nline;

%% Line parameters and adding extra line

% lines.TXnum = [0 lines.TXnum];
% lines.TXname = ['SUB' lines.TXname];
% 
% lines.RXnum = [1 lines.RXnum];
% lines.RXname = ['A1' lines.RXname];
% 
% lines.PH = [ones(3,1) lines.PH];
% 
% lines.config = ['12' lines.config];
% 
% lines.length = [3.048 lines.length];
% 
% lines.nline = lines.nline + 1;
% 
% lines.FZ(:,:,2:end+1) = lines.FZ;
% lines.FZ(:,:,1) = zeros(3,3);
% 
% lines.FZpu(:,:,2:end+1) = lines.FZpu;
% lines.FZpu(:,:,1) = zeros(3,3);
% 
% nline = lines.nline;

%% Load parameters

loads.spu = 2*loads.spu;

%% Controller parameters

controllers.wmaxpu = 1*controllers.wmaxpu

%% Simulation parameters

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim.Vnom = Vnom;

sim.Vfbs = [];
sim.Lfbs = [];
sim.Hfbs = [];

sim.rho = 0.1;

%%

% OPTSOL = Solver_LinDist3Flow_TX_minP0_20160603_mag(feeder, nodes, lines, configs, loads, caps, controllers, sim)
% 
% nvar = OPTSOL.nvar;
% Xa = OPTSOL.X(1:nvar);
% Xb = OPTSOL.X(nvar+1:2*nvar);
% Xc = OPTSOL.X(2*nvar+1:3*nvar);
% 
% XX = [Xa Xb Xc]';
% 
% Yopt = XX(:,1:nnode);
% XX(:,1:nnode) = [];
% 
% Yopt
% Vopt = sqrt(Yopt)
% 
% Popt = XX(:,1:nline);
% XX(:,1:nline) = [];
% 
% Qopt = XX(:,1:nline);
% XX(:,1:nline) = [];
% 
% Sopt = Popt + 1j*Qopt
% 
% uopt = XX(:,1:nnode);
% XX(:,1:nnode) = [];
% 
% vopt = XX(:,1:nnode);
% XX(:,1:nnode) = [];
% 
% wopt = uopt + 1j*vopt
% 
% dem = loads.spu.*(loads.aPQ + loads.aZ.*Yopt).*nodes.PH
% 
% for k1 = 2:nnode
%     sopt(:,k1) = sum(Sopt(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
%         sum(Sopt(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2);
% end
% sopt

%%

OPTSOL = Solver_LinDist3Flow_TX_minP0_20160603_magangle(feeder, nodes, lines, configs, loads, caps, controllers, sim)

nvar = OPTSOL.nvar;
Xa = OPTSOL.X(1:nvar);
Xb = OPTSOL.X(nvar+1:2*nvar);
Xc = OPTSOL.X(2*nvar+1:3*nvar);

XX = [Xa Xb Xc]';

Yopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

Yopt
Vopt = sqrt(Yopt)

Dopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

Dopt

Popt = XX(:,1:nline);
XX(:,1:nline) = [];

Qopt = XX(:,1:nline);
XX(:,1:nline) = [];

Sopt = Popt + 1j*Qopt

uopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

vopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

wopt = uopt + 1j*vopt

dem = loads.spu.*(loads.aPQ + loads.aZ.*Yopt).*nodes.PH

for k1 = 2:nnode
    sopt(:,k1) = sum(Sopt(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
        sum(Sopt(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2);
end
sopt


%%


% %% Obtain Dist3Flow solution
% 
% disp('----- Dist3Flow Formulation -----')
% 
% sim.iter = 0;
% sim.ramp = 0;
% optsol = Solver_LinDist3Flow_YBalance_20160603(feeder,loads,controllers,sim);
% 
% disp(['CVX Status: ' optsol.cvx_status])
% 
% nvar = optsol.nvar;
% Xa = optsol.X(1:nvar);
% Xb = optsol.X(nvar+1:2*nvar);
% Xc = optsol.X(2*nvar+1:3*nvar);
% 
% XX = [Xa Xb Xc]';
% 
% Yld3f = XX(:,1:n);
% Vmagapp = XX(:,1*n+1:2*n);
% Vmagld3f = sqrt(Yld3f);
% 
% Dld3f = XX(:,2*n+1:3*n);
% 
% Pld3f = XX(:,3*n+1:4*n);
% Qld3f = XX(:,4*n+1:5*n);
% Sld3f = Pld3f + j*Qld3f;
% 
% uld3f = XX(:,5*n+1:6*n);
% vld3f = XX(:,6*n+1:7*n);
% wld3f = uld3f + j*vld3f;
% 
% controllers.wpu = wld3f;
% FBScon = FBS_3phase_fun_20160603(feeder,loads,controllers,V0);
% 
% disp('/////////////////////////')
% disp(' ')
% 
% %% Plot feeder state with control
% 
% FBScon.V(feeder.PH == 0) = NaN;
% 
% % Plot voltage magnitude across feeder
% figure, box on, hold on
% for kpath = 1:size(feeder.paths,1)
%     temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
%     plot(temppath,abs(FBScon.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,abs(FBScon.V(2,temppath)),'gx','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,abs(FBScon.V(3,temppath)),'b.','MarkerSize',25,'LineWidth',2)
%     
% %     plot(temppath,abs(FBScon.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
% %     plot(temppath,abs(FBScon.V(2,temppath)),'g--','MarkerSize',12.5,'LineWidth',1.5)
% %     plot(temppath,abs(FBScon.V(3,temppath)),'b--','MarkerSize',25,'LineWidth',1.5)
% end
% plot(0.5:n+0.5,0.95*ones(n+1),'k--','MarkerSize',12,'LineWidth',2)
% set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'YTick',0.94:0.01:1.01,'FontSize',12,'FontWeight','bold')
% % title('Voltage Profile of Base (No Control) Case','FontSize',12,'FontWeight','bold')
% % legend({'a','b','c'},'FontSize',12,'FontWeight','bold','location','southwest')
% legend({'|V_{n}^{a}|','|V_{n}^{b}|','|V_{n}^{c}|'},'FontSize',12,'FontWeight','bold','location','southwest')
% xlabel('Node','FontSize',12,'FontWeight','bold')
% ylabel('Voltage Magnitude [pu]','FontSize',12,'FontWeight','bold')
% axis([0.5 n+0.5 0.94 1.00])
% pbaspect([1 0.3125 1])
% 
% print('-f3','-dpng','C:\Users\Michael\Desktop\temp\png\13node_balance_control.png')
% print('-f3','-depsc','C:\Users\Michael\Desktop\temp\eps\13node_balance_control.eps')
% 
% % Plot voltage angle across feeder
% figure, box on, hold on   
% for kpath = 1:size(feeder.paths,1)
%     temppath = feeder.paths(kpath,feeder.paths(kpath,:) ~=0);
%     plot(temppath,180/pi*angle(FBScon.V(1,temppath)),'r+','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,180/pi*angle(FBScon.V(2,temppath))+120,'gx','MarkerSize',12.5,'LineWidth',2)
%     plot(temppath,180/pi*angle(FBScon.V(3,temppath))-120,'b.','MarkerSize',25,'LineWidth',2)
%     
% %     plot(temppath,180/pi*angle(FBScon.V(1,temppath)),'r--','MarkerSize',12.5,'LineWidth',1.5)
% %     plot(temppath,180/pi*angle(FBScon.V(2,temppath))+120,'g--','MarkerSize',12.5,'LineWidth',1.5)
% %     plot(temppath,180/pi*angle(FBScon.V(3,temppath))-120,'b--','MarkerSize',25,'LineWidth',1.5)
% end
% set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'FontSize',12,'FontWeight','bold')
% % title('Feeder Voltage Angle of Base (No Control) Case','FontSize',12,'FontWeight','bold')
% legend({'a','b','c'},'FontSize',12,'FontWeight','bold','location','southwest')
% xlabel('Node','FontSize',12,'FontWeight','bold')
% ylabel('\theta_{n}^{\phi}','FontSize',12,'FontWeight','bold')
% axis([0.5 n+0.5 -2 2])
% pbaspect([1 0.3125 1])
% 
% %% Plot feeder voltage magnitude imbalance
% 
% % close all
% 
% % rows - a - b ; b - c; a - c
% % columns - node
% 
% Vdifbase = zeros(3,n);
% Vdifcon = zeros(3,n);
% Vdifbase_norm = zeros(n);
% Vdifcon_norm = zeros(n);
% for knode = 1:n
%     if feeder.PH(1,knode) == 1 && feeder.PH(2,knode) == 1 && feeder.PH(3,knode) == 0
%         Vdifbase(1,knode) = abs(FBSbase.V(1,knode)) - abs(FBSbase.V(2,knode));
%         Vdifcon(1,knode) = abs(FBScon.V(1,knode)) - abs(FBScon.V(2,knode));
%     elseif feeder.PH(1,knode) == 0 && feeder.PH(2,knode) == 1 && feeder.PH(3,knode) == 1
%         Vdifbase(2,knode) = abs(FBSbase.V(2,knode)) - abs(FBSbase.V(3,knode));
%         Vdifcon(2,knode) = abs(FBScon.V(2,knode)) - abs(FBScon.V(3,knode));
%     elseif feeder.PH(1,knode) == 1 && feeder.PH(2,knode) == 0 && feeder.PH(3,knode) == 1
%         Vdifbase(3,knode) = abs(FBSbase.V(1,knode)) - abs(FBSbase.V(3,knode));
%         Vdifcon(3,knode) = abs(FBScon.V(1,knode)) - abs(FBScon.V(3,knode));
%     elseif feeder.PH(1,knode) == 1 && feeder.PH(2,knode) == 1 && feeder.PH(3,knode) == 1
%         Vdifbase(1,knode) = abs(FBSbase.V(1,knode)) - abs(FBSbase.V(2,knode));
%         Vdifbase(2,knode) = abs(FBSbase.V(2,knode)) - abs(FBSbase.V(3,knode));
%         Vdifbase(3,knode) = abs(FBSbase.V(1,knode)) - abs(FBSbase.V(3,knode));
%         Vdifcon(1,knode) = abs(FBScon.V(1,knode)) - abs(FBScon.V(2,knode));
%         Vdifcon(2,knode) = abs(FBScon.V(2,knode)) - abs(FBScon.V(3,knode));
%         Vdifcon(3,knode) = abs(FBScon.V(1,knode)) - abs(FBScon.V(3,knode));
%     end
%     Vdifbase_2norm(knode) = norm(Vdifbase(:,knode),2);
%     Vdifcon_2norm(knode) = norm(Vdifcon(:,knode),2);
% end
% 
% % Vdifbase_2norm(Vdifbase_2norm == 0) = NaN;
% % Vdifcon_2norm(Vdifcon_2norm == 0) = NaN;
% 
% Jbase = sum(Vdifbase_2norm)
% Jcon = sum(Vdifcon_2norm)
% 
% figure, box on
% semilogy(1:n,Vdifbase_2norm,'r*','MarkerSize',12.5,'LineWidth',2), hold on
% semilogy(1:n,Vdifcon_2norm,'g.','MarkerSize',20,'LineWidth',2), hold on
% set(gca,'XTick',1:n,'XTickLabel',feeder.nodelist,'YTick',logspace(-4,-1,4),'FontSize',12,'FontWeight','bold')
% % title(['Voltage Magnitude Imbalance'],'FontSize',12,'FontWeight','bold')
% legend({'Base','Control'},'FontSize',12,'FontWeight','bold','location','east')
% xlabel('Node','FontSize',12,'FontWeight','bold')
% ylabel('J_{k}','FontSize',12,'FontWeight','bold')
% axis([0.5 n+0.5 1e-4 1e-1])
% pbaspect([1 0.3125 1])
% 
% print('-f5','-dpng','C:\Users\Michael\Desktop\temp\png\13node_unbalance.png')
% print('-f5','-depsc','C:\Users\Michael\Desktop\temp\eps\13node_unbalance.eps')
% 
% %% Plot control effort
% 
% % close all
% 
% figure, box on, hold on
% plot(1:length(controllers.cnodes),imag(wld3f(1,controllers.cnodes)),'r+','MarkerSize',12.5,'LineWidth',2)
% plot(1:length(controllers.cnodes),imag(wld3f(2,controllers.cnodes)),'gx','MarkerSize',12.5,'LineWidth',2)
% plot(1:length(controllers.cnodes),imag(wld3f(3,controllers.cnodes)),'b.','MarkerSize',25,'LineWidth',2)
% set(gca,'XTick',1:length(controllers.cnodes),'XTickLabel',feeder.nodelist(controllers.cnodes),'YTick',-0.08:0.02:0.08,'FontSize',12,'FontWeight','bold')
% % title(['DER Control Input'],'FontSize',12,'FontWeight','bold')
% legend({'v_{n}^{a}','v_{n}^{b}','v_{n}^{c}'},'FontSize',12,'FontWeight','bold','location','northwest')
% xlabel('Node','FontSize',12,'FontWeight','bold')
% ylabel('DER VAr dispatch [pu]','FontSize',12,'FontWeight','bold')
% axis([0 length(controllers.cnodes)+1 -0.08 0.08])
% pbaspect([1 0.3125 1])
% 
% print('-f6','-dpng','C:\Users\Michael\Desktop\temp\png\13node_der_dispatch.png')
% print('-f6','-depsc','C:\Users\Michael\Desktop\temp\eps\13node_der_dispatch.eps')
