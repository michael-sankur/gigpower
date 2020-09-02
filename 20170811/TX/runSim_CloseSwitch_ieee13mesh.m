% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

% if exist('FBS_3phase_fun_20160603','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
%     path(path,'C:\Users\Michael\Desktop\reconfiguration\20160603');
% end

% path(path,genpath('C:\Users\Michael\Desktop\reconfiguration\Mesh\GEN'));
% path(path,'C:\Users\Michael\Desktop\reconfiguration\Mesh\SDP');
% path(path,'C:\Users\Michael\Desktop\reconfiguration\Mesh\TX');

path(path,genpath('~/Dropbox/Unbalanced LinDistflow/20170811/'));

% if exist('mxlpsolve','file') == 0
%     path(path,'C:\Users\Michael\Desktop\mxlp');
% end
% 
% if exist('cvx_begin','file') == 0
%     cd C:\Users\Michael\Desktop\cvx-w64\cvx
%     cvx_setup
% end

%% Load feeder

name1 = 'ieee_13node_mesh_open';
% name1 = 'mesh_09node_switch_open';

% fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fp = '~/Dropbox/Unbalanced LinDistflow/20170811/Feeders/';
fn = [name1 '.txt'];

[feeder1, nodes1, lines1, configs1, loads1, caps1, controllers1] = feeder_mapper_TX_function_20170215(name1, fp, fn);

%% Feeder paramaters

radflag1 = 1;
for k1 = 1:size(nodes1.FM,1)
    if sum(nodes1.FM(k1,:) == -1) > 1
        radflag1 = 0;
    end    
end
radflag1

nnode1 = nodes1.nnode;
nline1 = lines1.nline;

%% Load parameters

loads1.aPQ = 0.85*ones(3,nnode1).*nodes1.PH;
loads1.aI = zeros(3,nnode1);
loads1.aZ = 0.15*ones(3,nnode1).*nodes1.PH;

loads1.spu(:,3:15) = 0.75*loads1.spu(:,3:15);
loads1.spu(:,16:28) = 1.5*loads1.spu(:,16:28);

%% Capacitor parameters

caps1.cappu = 1*caps1.cappu;
% caps1.cappu(:,3:15) = 0.5*caps1.cappu(:,3:15);
% caps1.cappu(:,16:28) = 0.5*caps1.cappu(:,16:28);

%% Controller parameters

controllers1.wmaxpu = 0.5*controllers1.wmaxpu;

%% Simulation parameters

Vnom1 = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim1.Vnom = Vnom1;

sim1.Vfbs = [];
sim1.Lfbs = [];
sim1.Hfbs = [];

sim1.rho = 0.1;

tn1 = 12;
tn2 = 25;

% tn1 = 6;
% tn2 = 7;

%%

[nvar1, Aineq1, bineq1, Aeq1, beq1] = createCVXMatrixCloseSwitch(feeder1,nodes1,lines1,configs1,loads1,caps1,controllers1,sim1);

cvx_begin quiet
    expressions Zmag1 Zang1 Zw1;
    variable X1(3*nvar1);
    Zmag1 = (X1(tn1) - X1(tn2))^2 + (X1(nvar1+tn1) - X1(nvar1+tn2))^2 + (X1(2*nvar1+tn1) - X1(2*nvar1+tn2))^2;
    Zang1 = (X1(nnode1+tn1) - X1(nnode1+tn2))^2 + (X1(nvar1+nnode1+tn1) - X1(nvar1+nnode1+tn2))^2 ...
        + (X1(2*nvar1+nnode1+tn1) - X1(2*nvar1+nnode1+tn2))^2;
    for ph = 1:3
        for k1 = 1:nnode1
            if controllers1.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar1 + nnode1 + nnode1 +  nline1 + nline1 + k1;
                kv = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + nnode1 + k1;
                Zw1 = Zw1 + X1(ku)^2 + X1(kv)^2;
            end
        end
    end
    minimize(1000*Zmag1 + 0*Zang1 + 1*Zw1)
    subject to;
    Aeq1 * X1 == beq1;
    Aineq1 * X1 <= bineq1;
%     for ph = 1:3
%         for k1 = 3:nline
%             kP = (ph-1)*nvar + nnode + nnode + k1;
%             kQ = (ph-1)*nvar + nnode + nnode + nline + k1;
%             norm(X([kP kQ]),2) <= 0.2
%         end
%     end
    for ph = 1:3
        for k1 = 1:nnode1
            ku = (ph-1)*nvar1 + nnode1 + nnode1 +  nline1 + nline1 + k1;
            kv = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + nnode1 + k1;
            norm(X1([ku kv]),2) <= 1.0*controllers1.wmaxpu(ph,k1)
        end
    end
cvx_end

[Eopt1, Dopt1, Vopt1, Iopt1, Sopt1, wopt1, demopt1, sopt1] = ...
    parse_CVX_output(X1,nvar1,feeder1,nodes1,lines1,configs1,loads1,caps1,controllers1);

%%

[nvar2, Aineq2, bineq2, Aeq2, beq2] = createCVXMatrixCloseSwitch(feeder1,nodes1,lines1,configs1,loads1,caps1,controllers1,sim1);

cvx_begin quiet
    expressions Zmag2 Zang2 Zw2;
    variable X2(3*nvar1);
    Zmag2 = (X2(tn1) - X2(tn2))^2 + (X2(nvar1+tn1) - X2(nvar1+tn2))^2 + (X2(2*nvar1+tn1) - X2(2*nvar1+tn2))^2;
    Zang2 = (X2(nnode1+tn1) - X2(nnode1+tn2))^2 + (X2(nvar1+nnode1+tn1) - X2(nvar1+nnode1+tn2))^2 ...
        + (X2(2*nvar1+nnode1+tn1) - X2(2*nvar1+nnode1+tn2))^2;
    for ph = 1:3
        for k1 = 1:nnode1
            if controllers1.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar1 + nnode1 + nnode1 +  nline1 + nline1 + k1;
                kv = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + nnode1 + k1;
                Zw2 = Zw2 + X2(ku)^2 + X2(kv)^2;
            end
        end
    end
    minimize(1000*Zmag2 + 1000*Zang2 + 1*Zw2)
    subject to;
    Aeq1 * X2 == beq1;
    Aineq1 * X2 <= bineq1;
%     for ph = 1:3
%         for k1 = 3:nline
%             kP = (ph-1)*nvar + nnode + nnode + k1;
%             kQ = (ph-1)*nvar + nnode + nnode + nline + k1;
%             norm(X([kP kQ]),2) <= 0.2
%         end
%     end
    for ph = 1:3
        for k1 = 1:nnode1
            ku = (ph-1)*nvar1 + nnode1 + nnode1 +  nline1 + nline1 + k1;
            kv = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + nnode1 + k1;
            norm(X2([ku kv]),2) <= 1.0*controllers1.wmaxpu(ph,k1)
        end
    end
cvx_end

[Eopt2, Dopt2, Vopt2, Iopt2, Sopt2, wopt2, demopt2, sopt2] = ...
    parse_CVX_output(X2,nvar1,feeder1,nodes1,lines1,configs1,loads1,caps1,controllers1);

%%

controllers1.wpu = zeros(3,nnode1);

[VNR0, INR0, STXNR0, SRXNR0] = NR3(feeder1,nodes1,lines1,configs1,loads1,caps1,controllers1,[],[],2,Vnom1);

controllers1.wpu = 1*wopt1;

[VNR1, INR1, STXNR1, SRXNR1] = NR3(feeder1,nodes1,lines1,configs1,loads1,caps1,controllers1,Vopt1,Iopt1,2,Vnom1);

controllers1.wpu = 1*wopt2;

[VNR2, INR2, STXNR2, SRXNR2] = NR3(feeder1,nodes1,lines1,configs1,loads1,caps1,controllers1,Vopt1,Iopt1,2,Vnom1);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load feeder

name2 = 'ieee_13node_mesh_closed';
% name2 = 'mesh_09node_switch_closed';

% fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fp = '~/Dropbox/Unbalanced LinDistflow/20170811/Feeders/';
fn = [name2 '.txt'];

[feeder2, nodes2, lines2, configs2, loads2, caps2, controllers2] = feeder_mapper_TX_function_20170215(name2, fp, fn);

%% Feeder paramaters

radflag1 = 1;
for k1 = 1:size(nodes2.FM,1)
    if sum(nodes2.FM(k1,:) == -1) > 1
        radflag1 = 0;
    end    
end
radflag1

nnode1 = nodes2.nnode;
nline1 = lines2.nline;

%% Load parameters

loads2.aPQ = 0.85*ones(3,nnode1).*nodes2.PH;
loads2.aI = zeros(3,nnode1);
loads2.aZ = 0.15*ones(3,nnode1).*nodes2.PH;

loads2.spu(:,3:15) = 0.5*loads2.spu(:,3:15);
loads2.spu(:,16:28) = 3*loads2.spu(:,16:28);

%% Capacitor parameters

caps2.cappu = 1*caps2.cappu;
caps2.cappu(:,3:15) = 1*caps2.cappu(:,3:15);
caps2.cappu(:,16:28) = 0.25*caps2.cappu(:,16:28);

%% Controller parameters

controllers2.wmaxpu = controllers2.wmaxpu;

%% Simulation parameters

Vnom2 = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim2.Vnom = Vnom1;

sim2.Vfbs = [];
sim2.Lfbs = [];
sim2.Hfbs = [];

sim2.rho = 0.1;

%%

controllers2.wpu = zeros(3,nnode1);

[VNR3, INR3, STXNR3, SRXNR3] = NR3(feeder2,nodes2,lines2,configs2,loads2,caps2,controllers2,[],[],2,Vnom2);

controllers2.wpu = wopt1;

[VNR4, INR4, STXNR4, SRXNR4] = NR3(feeder2,nodes2,lines2,configs2,loads2,caps2,controllers2,[],[],2,Vnom2);

controllers2.wpu = wopt2;

[VNR5, INR5, STXNR5, SRXNR5] = NR3(feeder2,nodes2,lines2,configs2,loads2,caps2,controllers2,[],[],2,Vnom2);

%%

clc

disp 'Open, No Control'
abs(abs(VNR0(:,tn1)) - abs(VNR0(:,tn2)))
180/pi*(angle(VNR0(:,tn1)) - angle(VNR0(:,tn2)))
SNR0 = VNR0(:,tn1).*conj(lines2.FYpu(:,:,end)*(VNR0(:,tn1) - VNR0(:,tn2)))

disp 'Open, Mag Only Control'
abs(abs(VNR1(:,tn1)) - abs(VNR1(:,tn2)))
180/pi*(angle(VNR1(:,tn1)) - angle(VNR1(:,tn2)))
SNR1 = VNR1(:,tn1).*conj(lines2.FYpu(:,:,end)*(VNR1(:,tn1) - VNR1(:,tn2)))

disp 'Open, Phasor Control'
abs(abs(VNR2(:,tn1)) - abs(VNR2(:,tn2)))
180/pi*(angle(VNR2(:,tn1)) - angle(VNR2(:,tn2)))
SNR2 = VNR2(:,tn1).*conj(lines2.FYpu(:,:,end)*(VNR2(:,tn1) - VNR2(:,tn2)))

disp 'Closed, No Control'
VNR3(:,tn1).*conj(lines2.FYpu(:,:,end)*(VNR3(:,tn1) - VNR3(:,tn2)))

disp 'Closed, Mag Only Control'
VNR4(:,tn1).*conj(lines2.FYpu(:,:,end)*(VNR4(:,tn1) - VNR4(:,tn2)))

disp 'Closed, Phasor Control'
VNR5(:,tn1).*conj(lines2.FYpu(:,:,end)*(VNR5(:,tn1) - VNR5(:,tn2)))

% Vtn1strA = sprintf('& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\',VNR0(1,tn1),VNR0(1,tn1)/1j,VNR1(1,tn1),VNR1(1,tn1)/1j,VNR2(1,tn1),VNR2(1,tn1)/1j)
% Vtn1strB = sprintf('& & $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\',VNR0(2,tn1),VNR0(2,tn1)/1j,VNR1(2,tn1),VNR1(2,tn1)/1j,VNR2(2,tn1),VNR2(2,tn1)/1j)
% Vtn1strC = sprintf('& & $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\',VNR0(3,tn1),VNR0(3,tn1)/1j,VNR1(3,tn1),VNR1(3,tn1)/1j,VNR2(3,tn1),VNR2(3,tn1)/1j)
% 
% Vtn2strA = sprintf('& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\',VNR0(1,tn2),VNR0(1,tn2)/1j,VNR1(1,tn2),VNR1(1,tn2)/1j,VNR2(1,tn2),VNR2(1,tn2)/1j)
% Vtn2strB = sprintf('& & $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\',VNR0(2,tn2),VNR0(2,tn2)/1j,VNR1(2,tn2),VNR1(2,tn2)/1j,VNR2(2,tn2),VNR2(2,tn2)/1j)
% Vtn2strC = sprintf('& & $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\',VNR0(3,tn2),VNR0(3,tn2)/1j,VNR1(3,tn2),VNR1(3,tn2)/1j,VNR2(3,tn2),VNR2(3,tn2)/1j)
% 
% magstrA = sprintf('& $a$ & %0.4f & %0.4f & %0.4f \\\\',(abs(VNR0(1,tn1)) - abs(VNR0(1,tn2))),(abs(VNR1(1,tn1)) - abs(VNR1(1,tn2))),(abs(VNR2(1,tn1)) - abs(VNR2(1,tn2))))
% magstrB = sprintf('& & $b$ & %0.4f & %0.4f & %0.4f \\\\',(abs(VNR0(2,tn1)) - abs(VNR0(2,tn2))),(abs(VNR1(2,tn1)) - abs(VNR1(2,tn2))),(abs(VNR2(2,tn1)) - abs(VNR2(2,tn2))))
% magstrC = sprintf('& & $c$ & %0.4f & %0.4f & %0.4f \\\\',(abs(VNR0(3,tn1)) - abs(VNR0(3,tn2))),(abs(VNR1(3,tn1)) - abs(VNR1(3,tn2))),(abs(VNR2(3,tn1)) - abs(VNR2(3,tn2))))
% 
% angstrA = sprintf('& $a$ & %0.4f & %0.4f & %0.4f \\\\',180/pi*(angle(VNR0(1,tn1)) - angle(VNR0(1,tn2))),180/pi*(angle(VNR1(1,tn1)) - angle(VNR1(1,tn2))),180/pi*(angle(VNR2(1,tn1)) - angle(VNR2(1,tn2))))
% angstrB = sprintf('& & $b$ & %0.4f & %0.4f & %0.4f \\\\',180/pi*(angle(VNR0(2,tn1)) - angle(VNR0(2,tn2))),180/pi*(angle(VNR1(2,tn1)) - angle(VNR1(2,tn2))),180/pi*(angle(VNR2(2,tn1)) - angle(VNR2(2,tn2))))
% angstrC = sprintf('& & $c$ & %0.4f & %0.4f & %0.4f \\\\',180/pi*(angle(VNR0(3,tn1)) - angle(VNR0(3,tn2))),180/pi*(angle(VNR1(3,tn1)) - angle(VNR1(3,tn2))),180/pi*(angle(VNR2(3,tn1)) - angle(VNR2(3,tn2))))
% 
% powstrA = sprintf('& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\',SNR0(1),SNR0(1)/1j,SNR1(1),SNR1(1)/1j,SNR2(1),SNR2(1)/1j)
% powstrB = sprintf('& & $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\',SNR0(2),SNR0(2)/1j,SNR1(2),SNR1(2)/1j,SNR2(2),SNR2(2)/1j)
% powstrC = sprintf('& & $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\',SNR0(3),SNR0(3)/1j,SNR1(3),SNR1(3)/1j,SNR2(3),SNR2(3)/1j)

%%

fid = fopen(['C:\Users\Michael\Desktop\temp\output\switch\' name1 '.txt'],'w');

Vtn1strA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(1,tn1),VNR0(1,tn1)/1j,VNR1(1,tn1),VNR1(1,tn1)/1j,VNR2(1,tn1),VNR2(1,tn1)/1j);
Vtn1strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(2,tn1),VNR0(2,tn1)/1j,VNR1(2,tn1),VNR1(2,tn1)/1j,VNR2(2,tn1),VNR2(2,tn1)/1j);
Vtn1strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(3,tn1),VNR0(3,tn1)/1j,VNR1(3,tn1),VNR1(3,tn1)/1j,VNR2(3,tn1),VNR2(3,tn1)/1j);
fprintf(fid,'\n');

Vtn1strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(1,tn1)),180/pi*angle(VNR0(1,tn1)),abs(VNR1(1,tn1)),180/pi*angle(VNR1(1,tn1)),abs(VNR2(1,tn1)),180/pi*angle(VNR2(1,tn1)));
Vtn1strA = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(2,tn1)),180/pi*angle(VNR0(2,tn1)),abs(VNR1(2,tn1)),180/pi*angle(VNR1(2,tn1)),abs(VNR2(2,tn1)),180/pi*angle(VNR2(2,tn1)));
Vtn1strA = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(3,tn1)),180/pi*angle(VNR0(3,tn1)),abs(VNR1(3,tn1)),180/pi*angle(VNR1(3,tn1)),abs(VNR2(3,tn1)),180/pi*angle(VNR2(3,tn1)));
fprintf(fid,'\n');

Vtn2strA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(1,tn2),VNR0(1,tn2)/1j,VNR1(1,tn2),VNR1(1,tn2)/1j,VNR2(1,tn2),VNR2(1,tn2)/1j);
Vtn2strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(2,tn2),VNR0(2,tn2)/1j,VNR1(2,tn2),VNR1(2,tn2)/1j,VNR2(2,tn2),VNR2(2,tn2)/1j);
Vtn2strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(3,tn2),VNR0(3,tn2)/1j,VNR1(3,tn2),VNR1(3,tn2)/1j,VNR2(3,tn2),VNR2(3,tn2)/1j);
fprintf(fid,'\n');

Vtn2strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(1,tn2)),180/pi*angle(VNR0(1,tn2)),abs(VNR1(1,tn2)),180/pi*angle(VNR1(1,tn2)),abs(VNR2(1,tn2)),180/pi*angle(VNR2(1,tn2)));
Vtn2strA = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(2,tn2)),180/pi*angle(VNR0(2,tn2)),abs(VNR1(2,tn2)),180/pi*angle(VNR1(2,tn2)),abs(VNR2(2,tn2)),180/pi*angle(VNR2(2,tn2)));
Vtn2strA = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(3,tn2)),180/pi*angle(VNR0(3,tn2)),abs(VNR1(3,tn2)),180/pi*angle(VNR1(3,tn2)),abs(VNR2(3,tn2)),180/pi*angle(VNR2(3,tn2)));
fprintf(fid,'\n');

magstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f & %0.4f \\\\\n',(abs(VNR0(1,tn1)) - abs(VNR0(1,tn2))),(abs(VNR1(1,tn1)) - abs(VNR1(1,tn2))),(abs(VNR2(1,tn1)) - abs(VNR2(1,tn2))));
magstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f & %0.4f \\\\\n',(abs(VNR0(2,tn1)) - abs(VNR0(2,tn2))),(abs(VNR1(2,tn1)) - abs(VNR1(2,tn2))),(abs(VNR2(2,tn1)) - abs(VNR2(2,tn2))));
magstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f & %0.4f \\\\\n',(abs(VNR0(3,tn1)) - abs(VNR0(3,tn2))),(abs(VNR1(3,tn1)) - abs(VNR1(3,tn2))),(abs(VNR2(3,tn1)) - abs(VNR2(3,tn2))));
fprintf(fid,'\n');

angstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f & %0.4f \\\\\n',180/pi*(angle(VNR0(1,tn1)) - angle(VNR0(1,tn2))),180/pi*(angle(VNR1(1,tn1)) - angle(VNR1(1,tn2))),180/pi*(angle(VNR2(1,tn1)) - angle(VNR2(1,tn2))));
angstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f & %0.4f \\\\\n',180/pi*(angle(VNR0(2,tn1)) - angle(VNR0(2,tn2))),180/pi*(angle(VNR1(2,tn1)) - angle(VNR1(2,tn2))),180/pi*(angle(VNR2(2,tn1)) - angle(VNR2(2,tn2))));
angstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f & %0.4f \\\\\n',180/pi*(angle(VNR0(3,tn1)) - angle(VNR0(3,tn2))),180/pi*(angle(VNR1(3,tn1)) - angle(VNR1(3,tn2))),180/pi*(angle(VNR2(3,tn1)) - angle(VNR2(3,tn2))));
fprintf(fid,'\n');

powstrA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',SNR0(1),SNR0(1)/1j,SNR1(1),SNR1(1)/1j,SNR2(1),SNR2(1)/1j);
powstrB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',SNR0(2),SNR0(2)/1j,SNR1(2),SNR1(2)/1j,SNR2(2),SNR2(2)/1j);
powstrC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',SNR0(3),SNR0(3)/1j,SNR1(3),SNR1(3)/1j,SNR2(3),SNR2(3)/1j);
fprintf(fid,'\n');

for k1 = 1:nodes1.nnode
    if sum(controllers1.wmaxpu(:,k1)) > 0

        fprintf(fid,char(nodes1.nodelist(k1)));
        for ph = 1:3
            if controllers1.wmaxpu(ph,k1) == 0
                fprintf(fid,' & 0');                
            else                
                fprintf(fid,' & %0.4f + j%0.4f',real(wopt2(ph,k1)),imag(wopt2(ph,k1)));
            end
        end
        fprintf(fid,' \\\\\n');
    end

end

fclose(fid);
