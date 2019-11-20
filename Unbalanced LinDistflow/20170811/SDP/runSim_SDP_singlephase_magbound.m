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

Vnom = [1;
    1*exp(j*240*pi/180);
    1*exp(j*120*pi/180)];
sim.Vnom = Vnom(1);

slacknode = 1;
sim.slacknode = 1;

sim.rho = 0.1;

sim.eps = 0.5;
sim.tn1 = 4;
sim.tn2 = 6;

sim.lb = -0.0005;
sim.ub = 0.0005;

%%

controllers.wpu = zeros(3,nnode);

[VNR0, INR0, STXNR0, SRXNR0] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,[],[],slacknode,Vnom);

[sdpsol,SDPMAT] = Solver_SDP_singlephase_magbound(feeder, nodes, lines, configs, loads, caps, controllers, sim)

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
lam = zeros(1,length(eps));
r1 = zeros(1,length(eps));

k1 = 1;
% for k1 = 1:length(vref)
    for k2 = 1:length(eps)
        
        sim.vref = vref(k1);
        sim.eps = eps(k2);
        
        [sdpsol,SDPMAT] = Solver_SDP_singlephase_magbound(feeder, nodes, lines, configs, loads, caps, controllers, sim);
        
        if contains(sdpsol.cvxstatus,'Solved') == 0
            feas(k1,k2) = 0;
        elseif contains(sdpsol.cvxstatus,'Solved') == 1
            feas(k1,k2) = 1;
        end

        [V, D] = eig(sdpsol.Xsdp);
        lam(k1,k2) = abs(D(6,6)/D(5,5));
        
    end
% end

r1 = lam >=1e6

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

%%

wopt = sdpsol.wsdp

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




