% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

clc, clear all, close all

% path(path,genpath('C:\Users\Michael\Desktop\reconfiguration\Mesh\GEN'));
path(path,genpath('/home/michael/Dropbox/Unbalanced LinDistflow/20180101/'));


% if exist('mxlpsolve','file') == 0
%     path(path,'C:\Users\Michael\Desktop\mxlp');
% end
% 
if exist('cvx_begin','file') == 0
%     cd C:\Users\Michael\Desktop\cvx-w64\cvx
%     cvx_setup
%     cd ~/cvx/
%     cvx_setup
end

%% Load feeder

networkname = '05node_fullphase_mesh';

fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180101/NETWORKS/';
fn = [networkname '.txt'];

[network1] = network_mapper_function_20180101(networkname, fp, fn);

%% Network paramaters

nnode = network1.nodes.nnode;
nline = network1.lines.nline;

%% Load parameters

network1.loads.aPQ = 0.85*ones(3,nnode).*network1.nodes.PH;
network1.loads.aI = zeros(3,nnode);
network1.loads.aZ = 0.15*ones(3,nnode).*network1.nodes.PH;

%% Capacitor parameters

network1.caps.cappu = 1*network1.caps.cappu;
% caps.cappu(:,3:15) = 0.5*caps.cappu(:,3:15);
% caps.cappu(:,16:28) = 0.5*caps.cappu(:,16:28);

%% Controller parameters

network1.cons.wmaxpu = 0.5*network1.cons.wmaxpu;

%% Simulation parameters

slackidx = 1;
sim.slackidx = slackidx;

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim.Vnom = Vnom;

sim.VNR = [];
sim.Lmn = [];
sim.Hmn = [];

rho = 0.1;
sim.rho = rho;

tn1 = 4;
tn2 = 6;

%% Run CVX

[nvar, Aineq, bineq, Aeq, beq] = create_CVX_Matrix_MagnitudeBalance_20180101(network1,sim);

cvx_begin quiet
    expressions Zbal Zw;
    variable Xbal(3*nvar);
    for k1 = 3:nnode
        if strcmp(network1.nodes.phases,'ab')
            Zbal = Zbal + (Xbal(k1) - Xbal(nvar + k1))^2;
        end
        if strcmp(network1.nodes.phases,'ac')
            Zbal = Zbal + (Xbal(k1) - Xbal(2*nvar + k1))^2;            
        end
        if strcmp(network1.nodes.phases,'bc')
            Zbal = Zbal + (Xbal(nvar + k1) - Xbal(2*nvar + k1))^2;  
        end
        if strcmp(network1.nodes.phases,'abc')
            Zbal = Zbal + (Xbal(k1) - Xbal(nvar + k1))^2 ...
                + (Xbal(k1) - Xbal(2*nvar + k1))^2 ...
                + (Xbal(nvar + k1) - Xbal(2*nvar + k1))^2; 
        end        
    end
%     for ph = 1:3
%         for k1 = 1:nnode
%             if cons.wmaxpu(ph,k1) > 0
%                 ku = (ph-1)*nvar + nnode + nnode +  nline + nline + k1;
%                 kv = (ph-1)*nvar + nnode + nnode + nline + nline + nnode + k1;
%                 Zw = Zw + Xbal(ku)^2 + Xbal(kv)^2;
%             end
%         end
%     end
    minimize(1*Zbal + rho*1*Zw)
    subject to;
    Aeq * Xbal == beq;
    Aineq * Xbal <= bineq;
%     for ph = 1:3
%         for k1 = 3:nline
%             kP = (ph-1)*nvar + nnode + nnode + k1;
%             kQ = (ph-1)*nvar + nnode + nnode + nline + k1;
%             norm(X([kP kQ]),2) <= 0.2
%         end
%     end
%     for ph = 1:3
%         for k1 = 1:nnode
%             ku = (ph-1)*nvar + nnode + nnode + nline + nline + k1;
%             kv = (ph-1)*nvar + nnode + nnode + nline + nline + nnode + k1;
%             norm(Xbal([ku kv]),2) <= 1.0*cons.wmaxpu(ph,k1)
%         end
%     end
cvx_end

[Eopt, Dopt, Vopt, Iopt, Sopt, wopt, demopt, sopt] = ...
    parse_CVX_output(Xbal,nvar,network1);


%%

network1.cons.wpu = zeros(3,nnode);

[VNR0, INR0, STXNR0, SRXNR0] = NR3(network1,[],[],slackidx,Vnom);

network1.cons.wpu = 1*wopt;

[VNR1, INR1, STXNR1, SRXNR1] = NR3(network1,Vopt,Iopt,slackidx,Vnom);

%%

BCIB = NaN*ones(1,nnode);
CCIB = NaN*ones(1,nnode);

for k1 = 2:nnode
    if strcmp(network1.nodes.phases,'ab')
        BCIB(k1) = abs(abs(VNR0(1,k1)) - abs(VNR0(2,k1)));
        CCIB(k1) = abs(abs(VNR1(1,k1)) - abs(VNR1(2,k1)));
    elseif strcmp(network1.nodes.phases,'ac')
        BCIB(k1) = abs(abs(VNR0(1,k1)) - abs(VNR0(3,k1)));
        CCIB(k1) = abs(abs(VNR1(1,k1)) - abs(VNR1(3,k1)));
    elseif strcmp(network1.nodes.phases,'bc')
        BCIB(k1) = abs(abs(VNR0(2,k1)) - abs(VNR0(3,k1)));
        CCIB(k1) = abs(abs(VNR1(2,k1)) - abs(VNR1(3,k1)));
    elseif strcmp(network1.nodes.phases,'abc')
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

% nodes.nodelist{1} = '\infty';

close all

VNR0PLOT = VNR0; VNR0PLOT(network1.nodes.PH == 0) = NaN;
VNR1PLOT = VNR1; VNR1PLOT(network1.nodes.PH == 0) = NaN;

figure, box on, hold on
plot(1:nnode,abs(VNR0PLOT(1,:)),'r+','MarkerSize',10,'LineWidth',2)
plot(1:nnode,abs(VNR0PLOT(2,:)),'gx','MarkerSize',10,'LineWidth',2)
plot(1:nnode,abs(VNR0PLOT(3,:)),'b.','MarkerSize',20,'LineWidth',2)
plot(0:nnode+1,0.95*ones(size(0:nnode+1)),'k--','LineWidth',2)
set(gca,'XTick',1:nnode,'XTickLabel',network1.nodes.nodelist,'YTick',0.94:0.01:1.01, ...
    'FontWeight','bold','FontSize',12)
legend({'a','b','c'},'FontWeight','bold','FontSize',12,'location','southwest')
title('Base Case - No Control','FontWeight','bold','FontSize',12)
xlabel('Node','FontWeight','bold','FontSize',12)
ylabel('Voltage Magnitude [p.u.]','FontWeight','bold','FontSize',12)
axis([0.5 nnode+0.5 0.94 1.01])
pbaspect([1 0.375 1])

% print('-f1','-depsc',['C:\Users\Michael\Desktop\temp\eps\balance\' feedername '_base.eps'])


figure, box on, hold on
plot(1:nnode,abs(VNR1PLOT(1,:)),'r+','MarkerSize',10,'LineWidth',2)
plot(1:nnode,abs(VNR1PLOT(2,:)),'gx','MarkerSize',10,'LineWidth',2)
plot(1:nnode,abs(VNR1PLOT(3,:)),'b.','MarkerSize',20,'LineWidth',2)
plot(0:nnode+1,0.95*ones(size(0:nnode+1)),'k--','LineWidth',2)
set(gca,'XTick',1:nnode,'XTickLabel',network1.nodes.nodelist,'YTick',0.94:0.01:1.01, ...
    'FontWeight','bold','FontSize',12)
legend({'a','b','c'},'FontWeight','bold','FontSize',12,'location','southwest')
title('Control Case','FontWeight','bold','FontSize',12)
xlabel('Node','FontWeight','bold','FontSize',12)
ylabel('Voltage Magnitude [p.u.]','FontWeight','bold','FontSize',12)
axis([0.5 nnode+0.5 0.94 1.01])
pbaspect([1 0.375 1])

% print('-f2','-depsc',['C:\Users\Michael\Desktop\temp\eps\balance\' feedername '_control.eps'])


figure, box on, hold on
plot(1:nnode,BCIB,'r+','MarkerSize',10,'LineWidth',2)
plot(1:nnode,CCIB,'gx','MarkerSize',10,'LineWidth',2)
set(gca,'XTick',1:nnode,'XTickLabel',network1.nodes.nodelist, ...
    'FontWeight','bold','FontSize',12)
legend({'Base Case','Control Case'},'FontWeight','bold','FontSize',12,'location','northwest')
title('Voltage Imbalance','FontWeight','bold','FontSize',12)
xlabel('Node','FontWeight','bold','FontSize',12)
ylabel('Voltage Imbalance [p.u.]','FontWeight','bold','FontSize',12)
axis([0.5 nnode+0.5 0 0.05])
pbaspect([1 0.375 1])

% print('-f3','-depsc',['C:\Users\Michael\Desktop\temp\eps\balance\' feedername '_imbalance.eps'])
