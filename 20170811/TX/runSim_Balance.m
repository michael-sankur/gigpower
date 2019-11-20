% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

% if exist('FBS_3phase_fun_20160603','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
%     path(path,'C:\Users\Michael\Desktop\reconfiguration\20160603');
% end

path(path,genpath('C:\Users\Michael\Desktop\reconfiguration\Mesh\GEN'));
path(path,'C:\Users\Michael\Desktop\reconfiguration\Mesh\SDP');
path(path,'C:\Users\Michael\Desktop\reconfiguration\Mesh\TX');

if exist('mxlpsolve','file') == 0
    path(path,'C:\Users\Michael\Desktop\mxlp');
end

if exist('cvx_begin','file') == 0
    cd C:\Users\Michael\Desktop\cvx-w64\cvx
    cvx_setup
end

%% Load feeder

feedername = 'ieee_13node_balance';
% feedername = 'ieee_37node_balance';
feedername = 'mesh_09node_balance';
% feedername = '05node_singlephase_radial';

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

%% Load parameters

loads.aPQ = 0.85*ones(3,nnode).*nodes.PH;
loads.aI = zeros(3,nnode);
loads.aZ = 0.15*ones(3,nnode).*nodes.PH;

%% Capacitor parameters

% caps.cappu = 1*caps.cappu;

%% Controller parameters

% controllers.wmaxpu = 1*controllers.wmaxpu;

%% IEEE 13 node params


% loads.spu = 1.125*loads.spu;
% 
% caps.cappu = 0*caps.cappu;
% 
% controllers.wmaxpu = 0.25*controllers.wmaxpu;
% 
% rho = 0.5;


%% Mesh 09 node params

loads.spu = 0.5*loads.spu;

caps.cappu = 1*caps.cappu;

controllers.wmaxpu = 0.5*controllers.wmaxpu;

rho = 0.5;

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

% rho = 0.5;
sim.rho = rho;

%%

[nvar, Aineq, bineq, Aeq, beq] = createCVXMatrixBalance(feeder,nodes,lines,configs,loads,caps,controllers,sim);

cvx_begin quiet
    expressions Zbal Zw;
    variable Xbal(3*nvar);
    for k1 = 3:nnode
        if sum(nodes.PH(:,k1) == [1 1 0]') == 3
            Zbal = Zbal + (Xbal(k1) - Xbal(nvar + k1))^2;
        end
        if sum(nodes.PH(:,k1) == [1 0 1]') == 3
            Zbal = Zbal + (Xbal(k1) - Xbal(2*nvar + k1))^2;            
        end
        if sum(nodes.PH(:,k1) == [0 1 1]') == 3
            Zbal = Zbal + (Xbal(nvar + k1) - Xbal(2*nvar + k1))^2;  
        end
        if sum(nodes.PH(:,k1) == [1 1 1]') == 3
            Zbal = Zbal + (Xbal(k1) - Xbal(nvar + k1))^2 ...
                + (Xbal(k1) - Xbal(2*nvar + k1))^2 ...
                + (Xbal(nvar + k1) - Xbal(2*nvar + k1))^2; 
        end        
    end
    for ph = 1:3
        for k1 = 1:nnode
            if controllers.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar + nnode + nnode +  nline + nline + k1;
                kv = (ph-1)*nvar + nnode + nnode + nline + nline + nnode + k1;
                Zw = Zw + Xbal(ku)^2 + Xbal(kv)^2;
            end
        end
    end
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
    for ph = 1:3
        for k1 = 1:nnode
            ku = (ph-1)*nvar + nnode + nnode + nline + nline + k1;
            kv = (ph-1)*nvar + nnode + nnode + nline + nline + nnode + k1;
            norm(Xbal([ku kv]),2) <= 1.0*controllers.wmaxpu(ph,k1)
        end
    end
cvx_end

[Eopt, Dopt, Vopt, Iopt, Sopt, wopt, demopt, sopt] = ...
    parse_CVX_output(Xbal,nvar,feeder,nodes,lines,configs,loads,caps,controllers);


%%

controllers.wpu = zeros(3,nnode);

[VNR0, INR0, STXNR0, SRXNR0] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,[],[],slacknode,Vnom);

controllers.wpu = 1*wopt;

[VNR1, INR1, STXNR1, SRXNR1] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,Vopt,Iopt,slacknode,Vnom);

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

nodes.nodelist{1} = '\infty';

close all

VNR0PLOT = VNR0; VNR0PLOT(nodes.PH == 0) = NaN;
VNR1PLOT = VNR1; VNR1PLOT(nodes.PH == 0) = NaN;

figure, box on, hold on
plot(1:nnode,abs(VNR0PLOT(1,:)),'r+','MarkerSize',10,'LineWidth',2)
plot(1:nnode,abs(VNR0PLOT(2,:)),'gx','MarkerSize',10,'LineWidth',2)
plot(1:nnode,abs(VNR0PLOT(3,:)),'b.','MarkerSize',20,'LineWidth',2)
plot(0:nnode+1,0.95*ones(size(0:nnode+1)),'k--','LineWidth',2)
set(gca,'XTick',1:nnode,'XTickLabel',nodes.nodelist,'YTick',0.94:0.01:1.01, ...
    'FontWeight','bold','FontSize',12)
legend({'a','b','c'},'FontWeight','bold','FontSize',12,'location','southwest')
title('Base Case - No Control','FontWeight','bold','FontSize',12)
xlabel('Node','FontWeight','bold','FontSize',12)
ylabel('Voltage Magnitude [p.u.]','FontWeight','bold','FontSize',12)
axis([0.5 nnode+0.5 0.94 1.01])
pbaspect([1 0.375 1])

print('-f1','-depsc',['C:\Users\Michael\Desktop\temp\eps\balance\' feedername '_base.eps'])


figure, box on, hold on
plot(1:nnode,abs(VNR1PLOT(1,:)),'r+','MarkerSize',10,'LineWidth',2)
plot(1:nnode,abs(VNR1PLOT(2,:)),'gx','MarkerSize',10,'LineWidth',2)
plot(1:nnode,abs(VNR1PLOT(3,:)),'b.','MarkerSize',20,'LineWidth',2)
plot(0:nnode+1,0.95*ones(size(0:nnode+1)),'k--','LineWidth',2)
set(gca,'XTick',1:nnode,'XTickLabel',nodes.nodelist,'YTick',0.94:0.01:1.01, ...
    'FontWeight','bold','FontSize',12)
legend({'a','b','c'},'FontWeight','bold','FontSize',12,'location','southwest')
title('Control Case','FontWeight','bold','FontSize',12)
xlabel('Node','FontWeight','bold','FontSize',12)
ylabel('Voltage Magnitude [p.u.]','FontWeight','bold','FontSize',12)
axis([0.5 nnode+0.5 0.94 1.01])
pbaspect([1 0.375 1])

print('-f2','-depsc',['C:\Users\Michael\Desktop\temp\eps\balance\' feedername '_control.eps'])


figure, box on, hold on
plot(1:nnode,BCIB,'r+','MarkerSize',10,'LineWidth',2)
plot(1:nnode,CCIB,'gx','MarkerSize',10,'LineWidth',2)
set(gca,'XTick',1:nnode,'XTickLabel',nodes.nodelist, ...
    'FontWeight','bold','FontSize',12)
legend({'Base Case','Control Case'},'FontWeight','bold','FontSize',12,'location','northwest')
title('Voltage Imbalance','FontWeight','bold','FontSize',12)
xlabel('Node','FontWeight','bold','FontSize',12)
ylabel('Voltage Imbalance [p.u.]','FontWeight','bold','FontSize',12)
axis([0.5 nnode+0.5 0 0.05])
pbaspect([1 0.375 1])

print('-f3','-depsc',['C:\Users\Michael\Desktop\temp\eps\balance\' feedername '_imbalance.eps'])

%%

nodes.nodelist{1} = '$\\infty$';

fid = fopen(['C:\Users\Michael\Desktop\temp\output\balance\' feedername '_output.txt'],'w');

for k1 = 1:nnode
    if sum(controllers.wmaxpu(:,k1)) > 0

        fprintf(fid,char(nodes.nodelist(k1)));
        for ph = 1:3
            if nodes.PH(ph,k1) == 0
                fprintf(fid,' & -');
            elseif nodes.PH(ph,k1) == 1
                if controllers.wmaxpu(ph,k1) == 0
                    fprintf(fid,' & 0');                
                else
                    if imag(wopt(ph,k1)) >= 0
%                         fprintf(fid,' & %0.4f + j%0.4f',100*real(wopt(ph,k1)),100*imag(wopt(ph,k1)));
                        fprintf(fid,' & %0.5f + j%0.5f',1*real(wopt(ph,k1)),1*imag(wopt(ph,k1)));
                    else
%                         fprintf(fid,' & %0.4f - j%0.4f',100*real(wopt(ph,k1)),-100*imag(wopt(ph,k1)));
                        fprintf(fid,' & %0.5f - j%0.5f',1*real(wopt(ph,k1)),-1*imag(wopt(ph,k1)));
                    end
                end
            end
        end
        fprintf(fid,' \\\\\n');
    end

end

fclose(fid);

%%

fid = fopen(['C:\Users\Michael\Desktop\temp\output\balance\' feedername '_param.txt'],'w');

fprintf(fid,'nodes.phases\n\n');

fprintf(fid,'Node & ');
for k1 = 1:nnode
    
    fprintf(fid,[char(nodes.nodelist(k1))]);
    if k1 == nnode
        fprintf(fid,' \\\\');
    else
        fprintf(fid,' & ');
    end
        
end
fprintf(fid,'\n');
fprintf(fid,'$\\mathcal{P}_{m}$ & ');
for k1 = 1:nnode
    fprintf(fid,'$\\left\\{ ');
    if nodes.PH(1,k1) == 1 && nodes.PH(2,k1) == 0 && nodes.PH(3,k1) == 0
        fprintf(fid,'a');
    elseif nodes.PH(1,k1) == 0 && nodes.PH(2,k1) == 1 && nodes.PH(3,k1) == 0
        fprintf(fid,'b');
    elseif nodes.PH(1,k1) == 0 && nodes.PH(2,k1) == 0 && nodes.PH(3,k1) == 1
        fprintf(fid,'c');
    elseif nodes.PH(1,k1) == 1 && nodes.PH(2,k1) == 1 && nodes.PH(3,k1) == 0
        fprintf(fid,'a,b');
    elseif nodes.PH(1,k1) == 1 && nodes.PH(2,k1) == 0 && nodes.PH(3,k1) == 1
        fprintf(fid,'a,c');
    elseif nodes.PH(1,k1) == 0 && nodes.PH(2,k1) == 1 && nodes.PH(3,k1) == 1
        fprintf(fid,'b,c');
    elseif nodes.PH(1,k1) == 1 && nodes.PH(2,k1) == 1 && nodes.PH(3,k1) == 1
        fprintf(fid,'a,b,c');       
    end    
    fprintf(fid,' \\right\\}$');
    if k1 == nnode
        fprintf(fid,' \\\\');
    else
        fprintf(fid,' & ');
    end
end

FZpu = 100*lines.FZpu;

fprintf(fid,'\n\nlines.phases\n\n');
for k1 = 1:nline
    
    fprintf(fid,char(nodes.nodelist(lines.TXnum(k1))));
    fprintf(fid,' & ');
    fprintf(fid,char(nodes.nodelist(lines.RXnum(k1))));
    fprintf(fid,' & ');
    fprintf(fid,'$\\left\\{ ');
    if lines.PH(1,k1) == 1 && lines.PH(2,k1) == 0 && lines.PH(3,k1) == 0
        fprintf(fid,'a');
    elseif lines.PH(1,k1) == 0 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 0
        fprintf(fid,'b');
    elseif lines.PH(1,k1) == 0 && lines.PH(2,k1) == 0 && lines.PH(3,k1) == 1
        fprintf(fid,'c');
    elseif lines.PH(1,k1) == 1 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 0
        fprintf(fid,'a,b');
    elseif lines.PH(1,k1) == 1 && lines.PH(2,k1) == 0 && lines.PH(3,k1) == 1
        fprintf(fid,'a,c');
    elseif lines.PH(1,k1) == 0 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 1
        fprintf(fid,'b,c');
    elseif lines.PH(1,k1) == 1 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 1
        fprintf(fid,'a,b,c');       
    end
    fprintf(fid,' \\right\\}$');
    fprintf(fid,' & ');
%     fprintf(fid,'$\\mathbf{Z}_{mn} = \\begin{bmatrix} ');
    fprintf(fid,'$\\begin{bmatrix} ');
    if lines.PH(1,k1) == 1 && lines.PH(2,k1) == 0 && lines.PH(3,k1) == 0
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(1,1,k1)),imag(FZpu(1,1,k1)));
    elseif lines.PH(1,k1) == 0 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 0
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(2,2,k1)),imag(FZpu(2,2,k1)));
    elseif lines.PH(1,k1) == 0 && lines.PH(2,k1) == 0 && lines.PH(3,k1) == 1
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(3,3,k1)),imag(FZpu(3,3,k1)));
    elseif lines.PH(1,k1) == 1 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 0
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(1,1,k1)),imag(FZpu(1,1,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(1,2,k1)),imag(FZpu(1,2,k1)));
        fprintf(fid,' \\\\ ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(2,1,k1)),imag(FZpu(2,1,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(2,2,k1)),imag(FZpu(2,2,k1)));
    elseif lines.PH(1,k1) == 1 && lines.PH(2,k1) == 0 && lines.PH(3,k1) == 1
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(1,1,k1)),imag(FZpu(1,1,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(1,3,k1)),imag(FZpu(1,3,k1)));
        fprintf(fid,' \\\\ ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(3,1,k1)),imag(FZpu(3,1,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(3,3,k1)),imag(FZpu(3,3,k1)));
    elseif lines.PH(1,k1) == 0 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 1
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(2,2,k1)),imag(FZpu(2,2,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(2,3,k1)),imag(FZpu(2,3,k1)));
        fprintf(fid,' \\\\ ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(3,2,k1)),imag(FZpu(3,2,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(3,3,k1)),imag(FZpu(3,3,k1)));
    elseif lines.PH(1,k1) == 1 && lines.PH(2,k1) == 1 && lines.PH(3,k1) == 1
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(1,1,k1)),imag(FZpu(1,1,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(1,2,k1)),imag(FZpu(1,2,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(1,3,k1)),imag(FZpu(1,3,k1)));
        fprintf(fid,' \\\\ ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(2,1,k1)),imag(FZpu(2,1,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(2,2,k1)),imag(FZpu(2,2,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(2,3,k1)),imag(FZpu(2,3,k1)));
        fprintf(fid,' \\\\ ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(3,1,k1)),imag(FZpu(3,1,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(3,2,k1)),imag(FZpu(3,2,k1)));
        fprintf(fid,' & ');
        fprintf(fid,'%0.4f + j %0.4f',real(FZpu(3,3,k1)),imag(FZpu(3,3,k1)));
    end
    fprintf(fid,' \\end{bmatrix}$');
    fprintf(fid,' \\\\\n');
    fprintf(fid,'\\hline\n');
    
end

fprintf(fid,'\n\nloads.spu\n\n');

for k1 = 1:nnode
    if sum(abs(loads.spu(:,k1))) > 0        
        fprintf(fid,char(nodes.nodelist(k1)));
        for ph = 1:3
            if nodes.PH(ph,k1) == 0
                fprintf(fid,' & -');
            elseif nodes.PH(ph,k1) == 1
                if loads.spu(ph,k1) == 0
                    fprintf(fid,' & 0');                
                else
                    if imag(loads.spu(ph,k1)) >= 0
                        fprintf(fid,' & %0.4f + j%0.4f',real(loads.spu(ph,k1)),imag(loads.spu(ph,k1)));
                    else
                        fprintf(fid,' & %0.4f - j%0.4f',real(loads.spu(ph,k1)),imag(loads.spu(ph,k1)));
                    end
                end
            end
        end
        fprintf(fid,' \\\\\n');        
    end
end

fprintf(fid,'\n\ncaps.cappu\n\n');


for k1 = 1:nnode
    if sum(abs(caps.cappu(:,k1))) > 0
        
        fprintf(fid,char(nodes.nodelist(k1)));
        for ph = 1:3
            if nodes.PH(ph,k1) == 0
                fprintf(fid,' & -');
            elseif nodes.PH(ph,k1) == 1
                if caps.cappu(ph,k1) == 0
                    fprintf(fid,' & 0');                
                else
                    fprintf(fid,' & j%0.4f',caps.cappu(ph,k1));
                end
            end
        end
        fprintf(fid,' \\\\\n');        
    end
end


fprintf(fid,'\n\ncontrollers.wmaxpu\n\n');


for k1 = 1:nnode
    if sum(abs(controllers.wmaxpu(:,k1))) > 0
        
        fprintf(fid,char(nodes.nodelist(k1)));
        for ph = 1:3
            if nodes.PH(ph,k1) == 0
                fprintf(fid,' & -');
            elseif nodes.PH(ph,k1) == 1
                if controllers.wmaxpu(ph,k1) == 0
                    fprintf(fid,' & 0');                
                else
                    fprintf(fid,' & %0.4f',controllers.wmaxpu(ph,k1));
                end
            end
        end
        fprintf(fid,' \\\\\n');        
    end
end

fclose(fid)
