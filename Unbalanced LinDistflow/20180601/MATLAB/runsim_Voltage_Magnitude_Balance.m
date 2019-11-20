% Michael Sankur - msankur@lbl.gov
% 2018.06.01

clc, clear all, close all

path(path,genpath('/home/michael/Dropbox/Unbalanced LinDistflow/20180601/MATLAB/'));

%% Load feeder

% Change the file path as necessary
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180601/NETWORKS/';
% Change file name as necessary
fn = '05node_fullphase_radial.txt';
fn = 'ieee_13node_balance.txt';

[network1] = network_mapper_function(fp, fn);

%% Network paramaters

nnode = network1.nodes.nnode;
nline = network1.lines.nline;

%% Load parameters

network1.loads.aPQ = 0.75*ones(3,nnode).*network1.nodes.PH;
network1.loads.aI = 0.10*ones(3,nnode).*network1.nodes.PH;
network1.loads.aZ = 0.15*ones(3,nnode).*network1.nodes.PH;

%% Capacitor parameters

network1.caps.cappu = 1*network1.caps.cappu;

%% Controller parameters

network1.cons.wmaxpu = 1*network1.cons.wmaxpu;

%% VVC parameters



%% Simulation parameters

slacknode = 1;
sim.slacknode = slacknode;

Vslack = 1*[1;
    1*exp(1j*-120*pi/180);
    1*exp(1j*120*pi/180)];
sim.Vslack = Vslack;

sim.VNR = [];
sim.Lmn = [];
sim.Hmn = [];

%% Run CVX

[nvar, Aineq, bineq, Aeq, beq] = create_linear_constraints_base(network1,slacknode,Vslack,sim);

cvx_begin quiet
    expressions Zbal Zw;
    variable Xbal(3*nvar);
    for k1 = 2:nnode
        if strcmp(network1.nodes.phases{k1},'ab')
            Zbal = Zbal + (Xbal(k1) - Xbal(nvar + k1))^2;
        end
        if strcmp(network1.nodes.phases{k1},'ac')
            Zbal = Zbal + (Xbal(k1) - Xbal(2*nvar + k1))^2;            
        end
        if strcmp(network1.nodes.phases{k1},'bc')
            Zbal = Zbal + (Xbal(nvar + k1) - Xbal(2*nvar + k1))^2;  
        end
        if strcmp(network1.nodes.phases{k1},'abc')
            Zbal = Zbal + (Xbal(k1) - Xbal(nvar + k1))^2 ...
                + (Xbal(k1) - Xbal(2*nvar + k1))^2 ...
                + (Xbal(nvar + k1) - Xbal(2*nvar + k1))^2; 
        end        
    end
    for ph = 1:3
        for k1 = 1:nnode
            if network1.cons.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar + nnode + nnode + nline + nline + k1;
                kv = (ph-1)*nvar + nnode + nnode + nline + nline + nnode + k1;
                Zw = Zw + Xbal(ku)^2 + Xbal(kv)^2;
            end
        end
    end
    minimize(1*Zbal + 0.1*Zw)
    subject to
    Aineq * Xbal <= bineq
    Aeq * Xbal == beq
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

[CVXopt, Eopt, Topt, Vopt, Iopt, Sopt, wopt, dopt, sopt] = ...
    parse_CVX_output(Xbal,nvar,network1);


%%

network1.cons.wpu = zeros(3,nnode);

[NRRES0, VNR0, INR0, STXNR0, SRXNR0] = NR3(network1,slacknode,Vslack,[],[],1e-9);

network1.cons.wpu = 1*wopt;

[NRRES1, VNR1, INR1, STXNR1, SRXNR1] = NR3(network1,slacknode,Vslack,Vopt,Iopt,1e-9);

%%

BCIB = NaN*ones(1,nnode);
CCIB = NaN*ones(1,nnode);

for k1 = 2:nnode
    if strcmp(network1.nodes.phases{k1},'ab')
        BCIB(k1) = abs(abs(VNR0(1,k1)) - abs(VNR0(2,k1)));
        CCIB(k1) = abs(abs(VNR1(1,k1)) - abs(VNR1(2,k1)));
    elseif strcmp(network1.nodes.phases{k1},'ac')
        BCIB(k1) = abs(abs(VNR0(1,k1)) - abs(VNR0(3,k1)));
        CCIB(k1) = abs(abs(VNR1(1,k1)) - abs(VNR1(3,k1)));
    elseif strcmp(network1.nodes.phases{k1},'bc')
        BCIB(k1) = abs(abs(VNR0(2,k1)) - abs(VNR0(3,k1)));
        CCIB(k1) = abs(abs(VNR1(2,k1)) - abs(VNR1(3,k1)));
    elseif strcmp(network1.nodes.phases{k1},'abc')
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
axis([0.5 nnode+0.5 0 0.1])
pbaspect([1 0.375 1])

% print('-f3','-depsc',['C:\Users\Michael\Desktop\temp\eps\balance\' feedername '_imbalance.eps'])
