% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

clc, clear all, close all

% path(path,genpath('C:\Users\Michael\Desktop\reconfiguration\Mesh\GEN'));
path(path,genpath('/home/michael/Dropbox/Unbalanced LinDistflow/20180101/'));


% if exist('mxlpsolve','file') == 0
%     path(path,'C:\Users\Michael\Desktop\mxlp');
% end
% 
% if exist('cvx_begin','file') == 0
%     cd C:\Users\Michael\Desktop\cvx-w64\cvx
%     cvx_setup
% end

%% Load feeder

networkname1 = 'ieee_13node_mesh_open';

fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180101/NETWORKS/';
fn = [networkname1 '.txt'];

[network1] = network_mapper_function_20180101(networkname1, fp, fn);

%% Network paramaters

nnode1 = network1.nodes.nnode;
nline1 = network1.lines.nline;

%% Load parameters

network1.loads.aPQ = 0.85*ones(3,nnode1).*network1.nodes.PH;
network1.loads.aI = zeros(3,nnode1);
network1.loads.aZ = 0.15*ones(3,nnode1).*network1.nodes.PH;

network1.loads.spu(:,3:15) = 0.75*network1.loads.spu(:,3:15);
network1.loads.spu(:,16:28) = 1.5*network1.loads.spu(:,16:28);

%% Capacitor parameters

network1.caps.cappu = 1*network1.caps.cappu;
% network1.caps.cappu(:,3:15) = 0.5*network1.caps.cappu(:,3:15);
% network1.caps.cappu(:,16:28) = 0.5*network1.caps.cappu(:,16:28);

%% Controller parameters

network1.cons.wmaxpu = 0.5*network1.cons.wmaxpu;

%% VVC parameters

% network1.vvc.qminpu = 0.5*network1.vvc.qminpu;
% network1.vvc.qmaxpu = 0.5*network1.vvc.qmaxpu;

%% Simulation parameters

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim.Vnom = Vnom;

sim.VNR = [];
sim.Lmn = [];
sim.Hmn = [];

sim.rho = 0.1;

tn1 = 12;
tn2 = 25;

%% Add line to close first switch

network2 = network1;
network2.lines.nline = network2.lines.nline+1;
network2.lines.TXnum(end+1) = tn1;
network2.lines.TXname(end+1) = network1.nodes.nodelist(tn1);
network2.lines.RXnum(end+1) = tn2;
network2.lines.RXname(end+1) = network1.nodes.nodelist(tn2);
lc = 3;
lm = 0.25;
network2.lines.phases{end+1} = network2.lines.phases{lc};
network2.lines.PH(:,end+1) = network2.lines.PH(:,lc);

network2.lines.config(end+1) = network2.lines.config(lc);
network2.lines.length(end+1) = lm*network2.lines.length(lc);

network2.lines.FZpu(:,:,end+1) = lm*network2.lines.FZpu(:,:,lc);
network2.lines.FRpu(:,:,end+1) = real(network2.lines.FZpu(:,:,end));
network2.lines.FXpu(:,:,end+1) = imag(network2.lines.FZpu(:,:,end));
network2.lines.FZ(:,:,end+1) = network2.lines.FZpu(:,:,end)*network2.base.Zbase;
network2.lines.FR(:,:,end+1) = network2.lines.FRpu(:,:,end)*network2.base.Zbase;
network2.lines.FX(:,:,end+1) = network2.lines.FXpu(:,:,end)*network2.base.Zbase;

network2.lines.FYpu(:,:,end+1) = 1/lm*network2.lines.FYpu(:,:,lc);
network2.lines.FGpu(:,:,end+1) = real(network2.lines.FZpu(:,:,end));
network2.lines.FBpu(:,:,end+1) = imag(network2.lines.FZpu(:,:,end));
network2.lines.FY(:,:,end+1) = network2.lines.FYpu(:,:,end)/network2.base.Zbase;
network2.lines.FG(:,:,end+1) = network2.lines.FGpu(:,:,end)/network2.base.Zbase;
network2.lines.FB(:,:,end+1) = network2.lines.FBpu(:,:,end)/network2.base.Zbase;

nnode2 = network2.nodes.nnode;
nline2 = network2.lines.nline;

%% Run CVX

[nvar1, Aineq1, bineq1, Aeq1, beq1] = create_CVX_Matrix_CloseSwitch_20180101(network1,sim);

cvx_begin quiet
    expressions Zmag1 Zang1 Zw1;
    variable X1(3*nvar1);
    Zmag1 = (X1(tn1) - X1(tn2))^2 + (X1(nvar1+tn1) - X1(nvar1+tn2))^2 + (X1(2*nvar1+tn1) - X1(2*nvar1+tn2))^2;
    Zang1 = (X1(nnode1+tn1) - X1(nnode1+tn2))^2 + (X1(nvar1+nnode1+tn1) - X1(nvar1+nnode1+tn2))^2 + (X1(2*nvar1+nnode1+tn1) - X1(2*nvar1+nnode1+tn2))^2;
    for ph = 1:3
        for k1 = 1:nnode1
            if network1.cons.wmaxpu(ph,k1) ~= 0
                ku = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + k1;
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
%     for ph = 1:3
%         for k1 = 1:nnode
%             ku = (ph-1)*nvar1 + nnode + nnode +  nline + nline + k1;
%             kv = (ph-1)*nvar1 + nnode + nnode + nline + nline + nnode + k1;
%             norm(X1([ku kv]),2) <= 1.0*cons.wmaxpu(ph,k1)
%         end
%     end
cvx_end

[Eopt1, Dopt1, Vopt1, Iopt1, Sopt1, wopt1, demopt1, sopt1] = parse_CVX_output(X1,nvar1,network1);

%% Run CVX 

[nvar2, Aineq2, bineq2, Aeq2, beq2] = create_CVX_Matrix_CloseSwitch_20180101(network1,sim);

cvx_begin quiet
    expressions Zmag2 Zang2 Zw2;
    variable X2(3*nvar1);
    Zmag2 = (X2(tn1) - X2(tn2))^2 + (X2(nvar2+tn1) - X2(nvar2+tn2))^2 + (X2(2*nvar2+tn1) - X2(2*nvar2+tn2))^2;
    Zang2 = (X2(nnode1+tn1) - X2(nnode1+tn2))^2 + (X2(nvar2+nnode1+tn1) - X2(nvar2+nnode1+tn2))^2 + (X2(2*nvar2+nnode1+tn1) - X2(2*nvar2+nnode1+tn2))^2;
    for ph = 1:3
        for k1 = 1:nnode1
            if network1.cons.wmaxpu(ph,k1) ~= 0
                ku = (ph-1)*nvar2 + nnode1 + nnode1 + nline1 + nline1 + k1;
                kv = (ph-1)*nvar2 + nnode1 + nnode1 + nline1 + nline1 + nnode1 + k1;
                Zw2 = Zw2 + X2(ku)^2 + X2(kv)^2;
            end
        end
    end
    minimize(1000*Zmag2 + 1000*Zang2 + 1*Zw2)
    subject to;
    Aeq2 * X2 == beq2;
    Aineq2 * X2 <= bineq2;
%     for ph = 1:3
%         for k1 = 3:nline
%             kP = (ph-1)*nvar + nnode + nnode + k1;
%             kQ = (ph-1)*nvar + nnode + nnode + nline + k1;
%             norm(X([kP kQ]),2) <= 0.2
%         end
%     end
%     for ph = 1:3
%         for k1 = 1:nnode
%             ku = (ph-1)*nvar2 + nnode + nnode +  nline + nline + k1;
%             kv = (ph-1)*nvar2 + nnode + nnode + nline + nline + nnode + k1;
%             norm(X2([ku kv]),2) <= 1.0*network1.cons.wmaxpu(ph,k1)
%         end
%     end
cvx_end

[Eopt2, Dopt2, Vopt2, Iopt2, Sopt2, wopt2, demopt2, sopt2] = parse_CVX_output(X2,nvar1,network1);

%%

network1.cons.wpu = zeros(3,nnode1);
network1.vvc.vvcpu = zeros(3,nnode1);
[VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, iter0] = NR3(network1,[],[],1,Vnom);
VNR0;
abs(VNR0);

network1.cons.wpu = wopt1;
network1.vvc.vvcpu = zeros(3,nnode1);
[VNR1, INR1, STXNR1, SRXNR1, iNR1, sNR1, iter1] = NR3(network1,Vopt1,Iopt1,1,Vnom);
VNR1;
abs(VNR1);

network1.cons.wpu = wopt2;
network1.vvc.vvcpu = zeros(3,nnode1);
[VNR2, INR2, STXNR2, SRXNR2, iNR2, sNR2, iter2] = NR3(network1,Vopt2,Iopt2,1,Vnom);
VNR2;
abs(VNR2);
VNR2(:,tn1) - VNR2(:,tn2)

%%

network1.vvc.vvcpu = zeros(3,nnode1);
for ph = 1:3
    for kn = 2:nnode1
        if network1.vvc.state(ph,kn) == 1
            qk = VVC_corrected(abs(VNR0(ph,kn)),network1.vvc.qminpu(ph,kn),network1.vvc.qmaxpu(ph,kn),network1.vvc.Vmin(ph,kn),network1.vvc.Vmax(ph,kn));
            network1.vvc.vvcpu(ph,kn) = qk;
        end        
    end
end
[VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, iter0] = NR3(network1,[],[],1,Vnom);
VNR0;
abs(VNR0);

network1.vvc.vvcpu = zeros(3,nnode1);
for ph = 1:3
    for kn = 2:nnode1
        if network1.vvc.state(ph,kn) == 1
            qk = VVC_corrected(abs(Vopt1(ph,kn)),network1.vvc.qminpu(ph,kn),network1.vvc.qmaxpu(ph,kn),network1.vvc.Vmin(ph,kn),network1.vvc.Vmax(ph,kn));
            network1.vvc.vvcpu(ph,kn) = qk;
        end        
    end
end
[VNR1, INR1, STXNR1, SRXNR1, iNR1, sNR1, iter1]= NR3(network1,[],[],1,Vnom);
VNR1;
abs(VNR1);

network1.vvc.vvcpu = zeros(3,nnode1);
for ph = 1:3
    for kn = 2:nnode1
        if network1.vvc.state(ph,kn) == 1
            qk = VVC_corrected(abs(Vopt2(ph,kn)),network1.vvc.qminpu(ph,kn),network1.vvc.qmaxpu(ph,kn),network1.vvc.Vmin(ph,kn),network1.vvc.Vmax(ph,kn));
            network1.vvc.vvcpu(ph,kn) = qk;
        end        
    end
end
[VNR2, INR2, STXNR2, SRXNR2, iNR2, sNR2, iter2]= NR3(network1,[],[],1,Vnom);
VNR2;
abs(VNR2);
VNR2(:,tn1) - VNR2(:,tn2)

%%

network1.cons.wpu = zeros(3,nnode1);

[VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, vvcNR0, iter0] = NR3VVC(network1,[],[],1,Vnom);
network1.vvc.vvcpu = vvcNR0;
VNR0;
abs(VNR0);

network1.cons.wpu = wopt1;
[VNR1, INR1, STXNR1, SRXNR1, iNR1, sNR1, vvcNR1, iter1] = NR3VVC(network1,Vopt1,Iopt1,1,Vnom);
network1.vvc.vvcpu = vvcNR1;
VNR1;
abs(VNR1);

network1.cons.wpu = wopt2;
[VNR2, INR2, STXNR2, SRXNR2, iNR2, sNR2, vvcNR2, iter2] = NR3VVC(network1,Vopt2,Iopt2,1,Vnom);
network1.vvc.vvcpu = vvcNR2;
VNR2;
abs(VNR2);
VNR2(:,tn1) - VNR2(:,tn2)


%%

disp 'Open, No Control'
VNR0(:,tn1) - VNR0(:,tn2)
abs(abs(VNR0(:,tn1)) - abs(VNR0(:,tn2)))
180/pi*(angle(VNR0(:,tn1)) - angle(VNR0(:,tn2)))
SNR0 = VNR0(:,tn1).*conj(network2.lines.FYpu(:,:,end)*(VNR0(:,tn1) - VNR0(:,tn2)))

disp 'Open, Mag Only Control'
VNR1(:,tn1) - VNR1(:,tn2)
abs(abs(VNR1(:,tn1)) - abs(VNR1(:,tn2)))
180/pi*(angle(VNR1(:,tn1)) - angle(VNR1(:,tn2)))
SNR1 = VNR1(:,tn1).*conj(network2.lines.FYpu(:,:,end)*(VNR1(:,tn1) - VNR1(:,tn2)))

disp 'Open, Phasor Control'
VNR2(:,tn1) - VNR2(:,tn2)
abs(abs(VNR2(:,tn1)) - abs(VNR2(:,tn2)))
180/pi*(angle(VNR2(:,tn1)) - angle(VNR2(:,tn2)))
SNR2 = VNR2(:,tn1).*conj(network2.lines.FYpu(:,:,end)*(VNR2(:,tn1) - VNR2(:,tn2)))

%%

% network2.cons.wpu = zeros(3,nnode1);
% 
% [VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, vvcNR0, iter0] = NR3VVC(network2,[],[],1,Vnom);
% network1.vvc.vvcpu = vvcNR0;
% VNR0;
% abs(VNR0);
% 
% network2.cons.wpu = wopt1;
% [VNR1, INR1, STXNR1, SRXNR1, iNR1, sNR1, vvcNR1, iter1] = NR3VVC(network2,Vopt1,[],1,Vnom);
% network1.vvc.vvcpu = vvcNR1;
% VNR1;
% abs(VNR1);
% 
% network2.cons.wpu = wopt2;
% [VNR2, INR2, STXNR2, SRXNR2, iNR2, sNR2, vvcNR2, iter2] = NR3VVC(network2,Vopt2,[],1,Vnom);
% network1.vvc.vvcpu = vvcNR2;
% VNR2;
% abs(VNR2);
% STXNR2(:,end)

%%

fid = fopen(['~/Desktop/temp/switch-' networkname1 '.txt'],'w');

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

for k1 = 1:network1.nodes.nnode
    if sum(network1.cons.wmaxpu(:,k1)) > 0

        fprintf(fid,char(network1.nodes.nodelist(k1)));
        for ph = 1:3
            if network1.cons.wmaxpu(ph,k1) == 0
                fprintf(fid,' & 0');                
            else                
                fprintf(fid,' & %0.4f + j%0.4f',real(wopt2(ph,k1)),imag(wopt2(ph,k1)));
            end
        end
        fprintf(fid,' \\\\\n');
    end

end

fclose(fid);
