% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

% if exist('FBS_3phase_fun_20160603','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
%     path(path,'C:\Users\Michael\Desktop\reconfiguration\20160603');
% end

path(path,'C:\Users\Michael\Desktop\reconfiguration\Mesh\GEN');
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
% feedername = '03node_fullphase';
feedername = '04node_fullphase';
feedername = '04node_multiphase';
feedername = '03node_fullphase';
feedername = '04node_singlephase';
feedername = '05node_singlephase_radial';
% feedername = '05node_singlephase_mesh';
% feedername = 'ieee_37node_singlephase';

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

%% Nodes parameters

%% Line parameters

%% Load parameters

% loads.aPQ = 0.85*ones(3,nnode).*nodes.PH;
% loads.aI = 0.0*ones(3,nnode).*nodes.PH;
% loads.aZ = 0.15*ones(3,nnode).*nodes.PH;

loads.spu = 1.5*loads.spu;

%% Capacitor parameters

caps.cappu = 1.0*caps.cappu;

%% Controller parameters

controllers.wmaxpu = 0.5*controllers.wmaxpu;

%% Simulation parameters

Vnom = 1;
sim.Vnom = Vnom;

slacknode = 1;
sim.slacknode = 1;

sim.rho = 0.1;

sim.eps = 0.0;
sim.tn = 3;
sim.vref = 0.980;

%%

Vnom = [1;
    1*exp(j*240*pi/180);
    1*exp(j*120*pi/180)];

controllers.wpu = zeros(3,nnode);

[VNR0, INR0, STXNR0, SRXNR0] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,[],[],slacknode,Vnom);

[sdpsol,SDPMAT] = Solver_SDP_singlephase_magref(feeder, nodes, lines, configs, loads, caps, controllers, sim)

[V, D] = eig(sdpsol.Xsdp);
lam = D(6,6)/D(5,5)

controllers.wpu = sdpsol.wsdp;
controllers.wpu(:,1) = 0;

[VNR1, INR1, STXNR1, SRXNR1] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,[],[],slacknode,Vnom);

abs(VNR1(1,:) - sdpsol.Xsdp(:,1).')

%%

vref = 0.95:0.005:1.00;
eps = 0:0.1:1.0;

feas = zeros(length(vref),length(eps));
lam = zeros(length(vref),length(eps));
r1 = zeros(length(vref),length(eps));

for k1 = 1:length(vref)
    for k2 = 1:length(eps)
        
        sim.vref = vref(k1);
        sim.eps = eps(k2);
        
        [sdpsol,SDPMAT] = Solver_SDP_singlephase_magref(feeder, nodes, lines, configs, loads, caps, controllers, sim);
        
        if contains(sdpsol.cvxstatus,'Solved') == 0
            feas(k1,k2) = 0;
        elseif contains(sdpsol.cvxstatus,'Solved') == 1
            feas(k1,k2) = 1;
        end

        [V, D] = eig(sdpsol.Xsdp);
        lam(k1,k2) = abs(D(6,6)/D(5,5));
        
%         if abs(D(6,6)/D(5,5)) >= 1e8
%             r1(k1,k2) = 1;
%         end
        
    end
end

%%

close all

r1 = lam >= 1e6;

figure, box on, hold on
for k1 = 1:length(vref)
    for k2 = 1:length(eps)
        
        if feas(k1,k2) == 0
            plot(vref(k1),eps(k2),'r.','MarkerSize',20)
        end
        
        if r1(k1,k2) == 0
            plot(vref(k1),eps(k2),'b.','MarkerSize',20)
        elseif r1(k1,k2) == 1
            plot(vref(k1),eps(k2),'g.','MarkerSize',20)
        end        
    end
end
set(gca,'XTick',0.95:0.01:1.0,'YTick',0:0.1:0.5,'FontWeight','bold','FontSize',12)
% legend({'A','B','C'},'FontWeight','bold','FontSize',12,'location','southwest')
title('Magnitude Reference - L1 Norm','FontWeight','bold','FontSize',12)
% title('Magnitude Reference - Quadratic','FontWeight','bold','FontSize',12)
xlabel('Voltage Reference \nu_{2} [p.u.]','FontWeight','bold','FontSize',12)
ylabel(char(949),'FontWeight','bold','FontSize',12)
axis([0.95 1.00 0 0.5])
pbaspect([1 0.375 1])

% print('-f1','-depsc',['C:\Users\Michael\Desktop\temp\eps\sdp\' feedername '_magref_L1_apppow.eps'])

print('-f1','-depsc',['C:\Users\Michael\Desktop\temp\eps\sdp\' feedername '_magref_Q2_apppow.eps'])

%% print results

fid = fopen(['C:\Users\Michael\Desktop\temp\output\sdp\' feedername '_magref.txt'],'w');

% fprintf(fid,'$V_{n}$ [p.u.] & ')
% for k1 = 1:nnode
%     fprintf(fid,'$%0.4f + j%0.4f$',real(VNR0(1,k1)),imag(VNR0(1,k1)));
%     if k1 == nnode
%         fprintf(fid,' \\\\');
%     else
%         fprintf(fid,' & ');
%     end    
% end
% 
% fprintf(fid,'\n')
% 
% fprintf(fid,'$V_{n}$ [p.u.] & ')
% for k1 = 1:nnode
%     fprintf(fid,'$%0.4f \\angle %0.4f^{\circ}$',abs(VNR0(1,k1)),180/pi*angle(VNR0(1,k1)));
%     if k1 == nnode
%         fprintf(fid,' \\\\');
%     else
%         fprintf(fid,' & ');
%     end    
% end

% fprintf(fid,'$V_{n}$ [p.u.] & ')
for k1 = 1:nnode
    fprintf(fid,nodes.nodelist{k1});
    fprintf(fid,' & ');
    fprintf(fid,'$%0.4f + j%0.4f$',real(VNR0(1,k1)),imag(VNR0(1,k1)));
    fprintf(fid,' & ');
    fprintf(fid,'$%0.4f \\angle %0.4f^{\\circ}$',abs(VNR0(1,k1)),180/pi*angle(VNR0(1,k1)));
    fprintf(fid,' \\\\\n');
end




