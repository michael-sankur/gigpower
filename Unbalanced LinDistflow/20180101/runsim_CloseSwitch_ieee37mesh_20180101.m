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

networkname = 'ieee_37node_mesh_open';

fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180101/NETWORKS/';
fn = [networkname '.txt'];

[network1] = network_mapper_function_20180101(networkname, fp, fn);

%% Network paramaters

nnode1 = network1.nodes.nnode;
nline1 = network1.lines.nline;

%% Load parameters

network1.loads.aPQ = 0.85*ones(3,nnode1).*network1.nodes.PH;
network1.loads.aI = zeros(3,nnode1);
network1.loads.aZ = 0.15*ones(3,nnode1).*network1.nodes.PH;

network1.loads.spu(:,3:39) = 1.5*network1.loads.spu(:,3:39);
network1.loads.spu(:,40:76) = 1.75*network1.loads.spu(:,40:76);

%% Capacitor parameters

network1.caps.cappu = 1*network1.caps.cappu;
% network1.caps.cappu(:,3:15) = 0.5*network1.caps.cappu(:,3:15);
% network1.caps.cappu(:,16:28) = 0.5*network1.caps.cappu(:,16:28);

%% Controller parameters

network1.cons.wmaxpu = 0.5*network1.cons.wmaxpu;

%% Simulation parameters

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim.Vnom = Vnom;

sim.VNR = [];
sim.Lmn = [];
sim.Hmn = [];

sim.rho = 0.1;

sw1 = [26 63];
sw2 = [15 52];

%% Add line to close first switch

network2 = network1;
network2.lines.nline = network2.lines.nline+1;
network2.lines.TXnum(end+1) = sw1(1);
network2.lines.TXname(end+1) = network1.nodes.nodelist(sw1(1));
network2.lines.RXnum(end+1) = sw1(2);
network2.lines.RXname(end+1) = network1.nodes.nodelist(sw1(2));
lc = 4;
lm = 4;
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

%% Add line to close second switch

network3 = network2;
network3.lines.nline = network3.lines.nline+1;
network3.lines.TXnum(end+1) = sw2(1);
network3.lines.TXname(end+1) = network1.nodes.nodelist(sw2(1));
network3.lines.RXnum(end+1) = sw2(2);
network3.lines.RXname(end+1) = network1.nodes.nodelist(sw2(2));
lc = 4;
lm = 4;
network3.lines.phases{end+1} = network3.lines.phases{lc};
network3.lines.PH(:,end+1) = network3.lines.PH(:,lc);

network3.lines.config(end+1) = network3.lines.config(lc);
network3.lines.length(end+1) = lm*network3.lines.length(lc);

network3.lines.FZpu(:,:,end+1) = lm*network3.lines.FZpu(:,:,lc);
network3.lines.FRpu(:,:,end+1) = real(network3.lines.FZpu(:,:,end));
network3.lines.FXpu(:,:,end+1) = imag(network3.lines.FZpu(:,:,end));
network3.lines.FZ(:,:,end+1) = network3.lines.FZpu(:,:,end)*network3.base.Zbase;
network3.lines.FR(:,:,end+1) = network3.lines.FRpu(:,:,end)*network3.base.Zbase;
network3.lines.FX(:,:,end+1) = network3.lines.FXpu(:,:,end)*network3.base.Zbase;

network3.lines.FYpu(:,:,end+1) = 1/lm*network3.lines.FYpu(:,:,lc);
network3.lines.FGpu(:,:,end+1) = real(network3.lines.FZpu(:,:,end));
network3.lines.FBpu(:,:,end+1) = imag(network3.lines.FZpu(:,:,end));
network3.lines.FY(:,:,end+1) = network3.lines.FYpu(:,:,end)/network3.base.Zbase;
network3.lines.FG(:,:,end+1) = network3.lines.FGpu(:,:,end)/network3.base.Zbase;
network3.lines.FB(:,:,end+1) = network3.lines.FBpu(:,:,end)/network3.base.Zbase;

nnode3 = network3.nodes.nnode;
nline3 = network3.lines.nline;

%% Simulate network1 with no control, both switches open

network1.cons.wpu = zeros(3,nnode1);
[VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, iter0] = NR3(network1,[],[],1,Vnom);
VD0 = VNR0(:,sw1(1)) - VNR0(:,sw1(2));
MD0 = abs(abs(VNR0(:,sw1(1))) - abs(VNR0(:,sw1(2))));
AD0 = 180/pi*(angle(VNR0(:,sw1(1))) - angle(VNR0(:,sw1(2))));
SS0 = VNR0(:,sw1(1)).*conj(network2.lines.FYpu(:,:,end)*(VNR0(:,sw1(1)) - VNR0(:,sw1(2))));

%% Run CVX for DER for closing first switch

[nvar1, Aineq1, bineq1, Aeq1, beq1] = create_CVX_Matrix_CloseSwitch_20180101(network1,sim);

cvx_begin quiet
    expressions Zmag1 Zang1 Zw1;
    variable X1(3*nvar1);
    Zmag1 = (X1(sw1(1)) - X1(sw1(2)))^2 + (X1(nvar1+sw1(1)) - X1(nvar1+sw1(2)))^2 + (X1(2*nvar1+sw1(1)) - X1(2*nvar1+sw1(2)))^2;
    Zang1 = (X1(nnode1+sw1(1)) - X1(nnode1+sw1(2)))^2 + (X1(nvar1+nnode1+sw1(1)) - X1(nvar1+nnode1+sw1(2)))^2 + (X1(2*nvar1+nnode1+sw1(1)) - X1(2*nvar1+nnode1+sw1(2)))^2;
    Zw1 = 0;
    for ph = 1:3
        for k1 = 1:nnode1
            if network1.cons.wmaxpu(ph,k1) ~= 0
                ku = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + k1;
                kv = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + nnode1 + k1;
                Zw1 = Zw1 + X1(ku)^2 + X1(kv)^2;
            end
        end
    end
    minimize(1000*Zmag1 + 1000*Zang1 + 1*Zw1)
    subject to
    Aeq1 * X1 == beq1;
    Aineq1 * X1 <= bineq1;
cvx_end

[Eopt1, Dopt1, Vopt1, Iopt1, Sopt1, wopt1, demopt1, sopt1] = parse_CVX_output(X1,nvar1,network1);

%% Simulation network1 with DER for closing first switch, both switches open

network1.cons.wpu = wopt1;
[VNR1, INR1, STXNR1, SRXNR1, iNR1, sNR1, iter1] = NR3(network1,Vopt1,Iopt1,1,Vnom);

VD1 = VNR1(:,sw1(1)) - VNR1(:,sw1(2));
MD1 = abs(abs(VNR1(:,sw1(1))) - abs(VNR1(:,sw1(2))));
AD1 = 180/pi*(angle(VNR1(:,sw1(1))) - angle(VNR1(:,sw1(2))));
SS1 = VNR1(:,sw1(1)).*conj(network2.lines.FYpu(:,:,end)*(VNR1(:,sw1(1)) - VNR1(:,sw1(2))));

%% Simulation network2 with DER for closing first switch, first switch closed

network2.cons.wpu = wopt1;
[VNR2, INR2, STXNR2, SRXNR2, iNR2, sNR2, iter2] = NR3(network2,Vopt1,[],1,Vnom);

VD2 = VNR2(:,sw1(1)) - VNR2(:,sw1(2));
MD2 = abs(abs(VNR2(:,sw1(1))) - abs(VNR2(:,sw1(2))));
AD2 = 180/pi*(angle(VNR2(:,sw1(1))) - angle(VNR2(:,sw1(2))));
SS2 = STXNR2(:,end);

%% Simulate network 2 with no control, first switch closed

network2.cons.wpu = zeros(3,nnode2);
[VNR3, INR3, STXNR3, SRXNR3, iNR3, sNR3, iter3] = NR3(network2,[],[],1,Vnom);

VD3 = VNR3(:,sw2(1)) - VNR3(:,sw2(2));
MD3 = abs(abs(VNR3(:,sw2(1))) - abs(VNR3(:,sw2(2))));
AD3 = 180/pi*(angle(VNR3(:,sw2(1))) - angle(VNR3(:,sw2(2))));
SS3 = VNR3(:,sw2(1)).*conj(network3.lines.FYpu(:,:,end)*(VNR3(:,sw2(1)) - VNR3(:,sw2(2))));

%% Run CVX for DER for closing second switch

[nvar2, Aineq2, bineq2, Aeq2, beq2] = create_CVX_Matrix_CloseSwitch_20180101(network2,sim);

cvx_begin quiet
    expressions Zmag2 Zang2 Zw2;
    variable X2(3*nvar2);
    Zmag2 = (X2(sw2(1)) - X2(sw2(2)))^2 + (X2(nvar2+sw2(1)) - X2(nvar2+sw2(2)))^2 + (X2(2*nvar2+sw2(1)) - X2(2*nvar2+sw2(2)))^2;
    Zang2 = (X2(nnode2+sw2(1)) - X2(nnode2+sw2(2)))^2 + (X2(nvar2+nnode2+sw2(1)) - X2(nvar2+nnode2+sw2(2)))^2 + (X2(2*nvar2+nnode2+sw2(1)) - X2(2*nvar2+nnode2+sw2(2)))^2;
    for ph = 1:3
        for k1 = 1:nnode2
            if network3.cons.wmaxpu(ph,k1) ~= 0
                ku = (ph-1)*nvar2 + nnode2 + nnode2 + nline2 + nline2 + k1;
                kv = (ph-1)*nvar2 + nnode2 + nnode2 + nline2 + nline2 + nnode2 + k1;
                Zw2 = Zw2 + X2(ku)^2 + X2(kv)^2;
            end
        end
    end
    minimize(1000*Zmag2 + 1000*Zang2 + 1*Zw2)
    subject to
    Aeq2 * X2 == beq2;
    Aineq2 * X2 <= bineq2;
cvx_end

[Eopt2, Dopt2, Vopt2, Iopt2, Sopt2, wopt2, demopt2, sopt2] = parse_CVX_output(X2,nvar2,network2);

%% Simulate network1 with DER for closing second switch, first switch closed

network2.cons.wpu = wopt2;
[VNR4, INR4, STXNR4, SRXNR4, iNR4, sNR4, iter4] = NR3(network2,Vopt2,Iopt2,1,Vnom);

VD4 = VNR4(:,sw2(1)) - VNR4(:,sw2(2));
MD4 = abs(abs(VNR4(:,sw2(1))) - abs(VNR4(:,sw2(2))));
AD4 = 180/pi*(angle(VNR4(:,sw2(1))) - angle(VNR4(:,sw2(2))));
SS4 = VNR4(:,sw2(1)).*conj(network3.lines.FYpu(:,:,end)*(VNR4(:,sw2(1)) - VNR4(:,sw2(2))));

%% Simulate network2 with DER for closing second switch, second switch closed

network3.cons.wpu = wopt2;
[VNR5, INR5, STXNR5, SRXNR5, iNR5, sNR5, iter5] = NR3(network3,[],[],1,Vnom);

VD5 = VNR5(:,sw2(1)) - VNR5(:,sw2(2));
MD5 = abs(abs(VNR5(:,sw2(1))) - abs(VNR5(:,sw2(2))));
AD5 = 180/pi*(angle(VNR5(:,sw2(1))) - angle(VNR5(:,sw2(2))));
SS5 = STXNR5(:,end);

%% Simulate network2 with no control, second switch closed

network3.cons.wpu = zeros(3,nnode2);
[VNR6, INR6, STXNR6, SRXNR6, iNR6, sNR6, iter6] = NR3(network3,[],[],1,Vnom);

VD6 = VNR6(:,sw2(1)) - VNR6(:,sw2(2));
MD6 = abs(abs(VNR6(:,sw2(1))) - abs(VNR6(:,sw2(2))));
AD6 = 180/pi*(angle(VNR6(:,sw2(1))) - angle(VNR6(:,sw2(2))));
SS6 = STXNR6(:,end);

%%

fid = fopen(['~/Desktop/temp/switch-' networkname '-sw1.txt'],'w');

Vtn1strA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(1,sw1(1)),VNR0(1,sw1(1))/1j,VNR1(1,sw1(1)),VNR1(1,sw1(1))/1j);
Vtn1strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(2,sw1(1)),VNR0(2,sw1(1))/1j,VNR1(2,sw1(1)),VNR1(2,sw1(1))/1j);
Vtn1strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(3,sw1(1)),VNR0(3,sw1(1))/1j,VNR1(3,sw1(1)),VNR1(3,sw1(1))/1j);
fprintf(fid,'\n');

Vtn1strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(1,sw1(1))),180/pi*angle(VNR0(1,sw1(1))),abs(VNR1(1,sw1(1))),180/pi*angle(VNR1(1,sw1(1))));
Vtn1strA = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(2,sw1(1))),180/pi*angle(VNR0(2,sw1(1))),abs(VNR1(2,sw1(1))),180/pi*angle(VNR1(2,sw1(1))));
Vtn1strA = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(3,sw1(1))),180/pi*angle(VNR0(3,sw1(1))),abs(VNR1(3,sw1(1))),180/pi*angle(VNR1(3,sw1(1))));
fprintf(fid,'\n');

Vtn1strA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(1,sw1(2)),VNR0(1,sw1(2))/1j,VNR1(1,sw1(2)),VNR1(1,sw1(2))/1j);
Vtn1strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(2,sw1(2)),VNR0(2,sw1(2))/1j,VNR1(2,sw1(2)),VNR1(2,sw1(2))/1j);
Vtn1strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(3,sw1(2)),VNR0(3,sw1(2))/1j,VNR1(3,sw1(2)),VNR1(3,sw1(2))/1j);
fprintf(fid,'\n');

Vtn2strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(1,sw1(2))),180/pi*angle(VNR0(1,sw1(2))),abs(VNR1(1,sw1(2))),180/pi*angle(VNR1(1,sw1(2))));
Vtn2strA = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(2,sw1(2))),180/pi*angle(VNR0(2,sw1(2))),abs(VNR1(2,sw1(2))),180/pi*angle(VNR1(2,sw1(2))));
Vtn2strA = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(3,sw1(2))),180/pi*angle(VNR0(3,sw1(2))),abs(VNR1(3,sw1(2))),180/pi*angle(VNR1(3,sw1(2))));
fprintf(fid,'\n');

magstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f \\\\\n',(abs(VNR0(1,sw1(1))) - abs(VNR0(1,sw1(2)))),(abs(VNR1(1,sw1(1))) - abs(VNR1(1,sw1(2)))));
magstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f \\\\\n',(abs(VNR0(2,sw1(1))) - abs(VNR0(2,sw1(2)))),(abs(VNR1(2,sw1(1))) - abs(VNR1(2,sw1(2)))));
magstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f \\\\\n',(abs(VNR0(3,sw1(1))) - abs(VNR0(3,sw1(2)))),(abs(VNR1(3,sw1(1))) - abs(VNR1(3,sw1(2)))));
fprintf(fid,'\n');

angstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f \\\\\n',180/pi*(angle(VNR0(1,sw1(1))) - angle(VNR0(1,sw1(2)))),180/pi*(angle(VNR1(1,sw1(1))) - angle(VNR1(1,sw1(2)))));
angstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f \\\\\n',180/pi*(angle(VNR0(2,sw1(1))) - angle(VNR0(2,sw1(2)))),180/pi*(angle(VNR1(2,sw1(1))) - angle(VNR1(2,sw1(2)))));
angstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f \\\\\n',180/pi*(angle(VNR0(3,sw1(1))) - angle(VNR0(3,sw1(2)))),180/pi*(angle(VNR1(3,sw1(1))) - angle(VNR1(3,sw1(2)))));
fprintf(fid,'\n');

powstrA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',SS0(1),SS0(1)/1j,SS1(1),SS1(1)/1j);
powstrB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',SS0(2),SS0(2)/1j,SS1(2),SS1(2)/1j);
powstrC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',SS0(3),SS0(3)/1j,SS1(3),SS1(3)/1j);
fprintf(fid,'\n');

for k1 = 1:network3.nodes.nnode
    if sum(network3.cons.wmaxpu(:,k1)) ~= 0

        fprintf(fid,char(network3.nodes.nodelist(k1)));
        for ph = 1:3
            if network3.cons.wmaxpu(ph,k1) == 0
                fprintf(fid,' & -');                
            else                
                fprintf(fid,' & %0.4f + j%0.4f',real(wopt1(ph,k1)),imag(wopt1(ph,k1)));
            end
        end
        fprintf(fid,' \\\\\n');
    end

end

fclose(fid);

%%

fid = fopen(['~/Desktop/temp/switch-' networkname '-sw2.txt'],'w');

Vtn1strA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    VNR3(1,sw2(1)),VNR3(1,sw2(1))/1j,VNR4(1,sw2(1)),VNR4(1,sw2(1))/1j);
Vtn1strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    VNR3(2,sw2(1)),VNR3(2,sw2(1))/1j,VNR4(2,sw2(1)),VNR4(2,sw2(1))/1j);
Vtn1strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    VNR3(3,sw2(1)),VNR3(3,sw2(1))/1j,VNR4(3,sw2(1)),VNR4(3,sw2(1))/1j);
fprintf(fid,'\n');

Vtn1strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR3(1,sw2(1))),180/pi*angle(VNR3(1,sw2(1))),abs(VNR4(1,sw2(1))),180/pi*angle(VNR4(1,sw2(1))));
Vtn1strA = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR3(2,sw2(1))),180/pi*angle(VNR3(2,sw2(1))),abs(VNR4(2,sw2(1))),180/pi*angle(VNR4(2,sw2(1))));
Vtn1strA = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR3(3,sw2(1))),180/pi*angle(VNR3(3,sw2(1))),abs(VNR4(3,sw2(1))),180/pi*angle(VNR4(3,sw2(1))));
fprintf(fid,'\n');

Vtn1strA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    VNR3(1,sw2(2)),VNR3(1,sw2(2))/1j,VNR4(1,sw2(2)),VNR4(1,sw2(2))/1j);
Vtn1strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    VNR3(2,sw2(2)),VNR3(2,sw2(2))/1j,VNR4(2,sw2(2)),VNR4(2,sw2(2))/1j);
Vtn1strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    VNR3(3,sw2(2)),VNR3(3,sw2(2))/1j,VNR4(3,sw2(2)),VNR4(3,sw2(2))/1j);
fprintf(fid,'\n');

Vtn2strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR3(1,sw2(2))),180/pi*angle(VNR3(1,sw2(2))),abs(VNR4(1,sw2(2))),180/pi*angle(VNR4(1,sw2(2))));
Vtn2strA = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR3(2,sw2(2))),180/pi*angle(VNR3(2,sw2(2))),abs(VNR4(2,sw2(2))),180/pi*angle(VNR4(2,sw2(2))));
Vtn2strA = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR3(3,sw2(2))),180/pi*angle(VNR3(3,sw2(2))),abs(VNR4(3,sw2(2))),180/pi*angle(VNR4(3,sw2(2))));
fprintf(fid,'\n');

magstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f \\\\\n',...
    (abs(VNR3(1,sw2(1))) - abs(VNR3(1,sw2(2)))),(abs(VNR4(1,sw2(1))) - abs(VNR4(1,sw2(2)))));
magstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f \\\\\n',...
    (abs(VNR3(2,sw2(1))) - abs(VNR3(2,sw2(2)))),(abs(VNR4(2,sw2(1))) - abs(VNR4(2,sw2(2)))));
magstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f \\\\\n',...
    (abs(VNR3(3,sw2(1))) - abs(VNR3(3,sw2(2)))),(abs(VNR4(3,sw2(1))) - abs(VNR4(3,sw2(2)))));
fprintf(fid,'\n');

angstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f \\\\\n',...
    180/pi*(angle(VNR3(1,sw2(1))) - angle(VNR3(1,sw2(2)))),180/pi*(angle(VNR4(1,sw2(1))) - angle(VNR4(1,sw2(2)))));
angstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f \\\\\n',...
    180/pi*(angle(VNR3(2,sw2(1))) - angle(VNR3(2,sw2(2)))),180/pi*(angle(VNR4(2,sw2(1))) - angle(VNR4(2,sw2(2)))));
angstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f \\\\\n',...
    180/pi*(angle(VNR3(3,sw2(1))) - angle(VNR3(3,sw2(2)))),180/pi*(angle(VNR4(3,sw2(1))) - angle(VNR4(3,sw2(2)))));
fprintf(fid,'\n');

powstrA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    SS3(1),SS3(1)/1j,SS4(1),SS4(1)/1j);
powstrB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    SS3(2),SS3(2)/1j,SS4(2),SS4(2)/1j);
powstrC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    SS3(3),SS3(3)/1j,SS4(3),SS4(3)/1j);
fprintf(fid,'\n');

for k1 = 1:network2.nodes.nnode
    if sum(network2.cons.wmaxpu(:,k1)) ~= 0
        
        fprintf(fid,char(network2.nodes.nodelist(k1)));
        for ph = 1:3
            if network2.cons.wmaxpu(ph,k1) == 0
                fprintf(fid,' & -');                
            else
                fprintf(fid,' & %0.4f + j%0.4f',real(wopt2(ph,k1)),imag(wopt2(ph,k1)));
            end% fid = fopen(['~/Desktop/temp/switch/' name1 '.txt'],'w');

        end
        fprintf(fid,' \\\\\n');
    end

end

fclose(fid);

%%

fid = fopen(['~/Desktop/temp/switch-' networkname '-sw1sw2.txt'],'w');

% Vtn1strA = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
%     VNR0(1,sw1(2)),VNR0(1,sw1(2))/1j,VNR1(1,sw1(2)),VNR1(1,sw1(2))/1j,VNR3(1,sw2(2)),VNR3(1,sw2(2))/1j,VNR4(1,sw2(2)),VNR4(1,sw2(2))/1j);
% Vtn1strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
%     VNR0(2,sw1(2)),VNR0(2,sw1(2))/1j,VNR1(2,sw1(2)),VNR1(2,sw1(2))/1j,VNR3(2,sw2(2)),VNR3(2,sw2(2))/1j,VNR4(2,sw2(2)),VNR4(2,sw2(2))/1j);
% Vtn1strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
%     VNR0(3,sw1(2)),VNR0(3,sw1(2))/1j,VNR1(3,sw1(2)),VNR1(3,sw1(2))/1j,VNR3(3,sw2(2)),VNR3(3,sw2(2))/1j,VNR4(3,sw2(2)),VNR4(3,sw2(2))/1j);
% fprintf(fid,'\n');

Vtn1strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(1,sw1(1))),180/pi*angle(VNR0(1,sw1(1))),abs(VNR1(1,sw1(1))),180/pi*angle(VNR1(1,sw1(1))),...
    abs(VNR3(1,sw2(1))),180/pi*angle(VNR3(1,sw2(1))),abs(VNR4(1,sw2(1))),180/pi*angle(VNR4(1,sw2(1))));
Vtn1strB = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(2,sw1(1))),180/pi*angle(VNR0(2,sw1(1))),abs(VNR1(2,sw1(1))),180/pi*angle(VNR1(2,sw1(1))),...
    abs(VNR3(1,sw2(1))),180/pi*angle(VNR3(2,sw2(1))),abs(VNR4(2,sw2(1))),180/pi*angle(VNR4(2,sw2(1))));
Vtn1strC = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(3,sw1(1))),180/pi*angle(VNR0(3,sw1(1))),abs(VNR1(3,sw1(1))),180/pi*angle(VNR1(3,sw1(1))),...
    abs(VNR3(1,sw2(1))),180/pi*angle(VNR3(3,sw2(1))),abs(VNR4(3,sw2(1))),180/pi*angle(VNR4(3,sw2(1))));
fprintf(fid,'\n');

% Vtn1strA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
%     VNR0(1,sw1(2)),VNR0(1,sw1(2))/1j,VNR1(1,sw1(2)),VNR1(1,sw1(2))/1j,VNR3(1,sw2(2)),VNR3(1,sw2(2))/1j,VNR4(1,sw2(2)),VNR4(1,sw2(2))/1j);
% Vtn1strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
%     VNR0(2,sw1(2)),VNR0(2,sw1(2))/1j,VNR1(2,sw1(2)),VNR1(2,sw1(2))/1j,VNR3(2,sw2(2)),VNR3(2,sw2(2))/1j,VNR4(2,sw2(2)),VNR4(2,sw2(2))/1j);
% Vtn1strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
%     VNR0(3,sw1(2)),VNR0(3,sw1(2))/1j,VNR1(3,sw1(2)),VNR1(3,sw1(2))/1j,VNR3(3,sw2(2)),VNR3(3,sw2(2))/1j,VNR4(3,sw2(2)),VNR4(3,sw2(2))/1j);
% fprintf(fid,'\n');

Vtn1strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(1,sw1(2))),180/pi*angle(VNR0(1,sw1(2))),abs(VNR1(1,sw1(2))),180/pi*angle(VNR1(1,sw1(2))),...
    abs(VNR3(1,sw2(2))),180/pi*angle(VNR3(1,sw2(2))),abs(VNR4(1,sw2(2))),180/pi*angle(VNR4(1,sw2(2))));
Vtn1strB = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(2,sw1(2))),180/pi*angle(VNR0(2,sw1(2))),abs(VNR1(2,sw1(2))),180/pi*angle(VNR1(2,sw1(2))),...
    abs(VNR3(1,sw2(2))),180/pi*angle(VNR3(2,sw2(2))),abs(VNR4(2,sw2(2))),180/pi*angle(VNR4(2,sw2(2))));
Vtn1strC = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
    abs(VNR0(3,sw1(2))),180/pi*angle(VNR0(3,sw1(2))),abs(VNR1(3,sw1(2))),180/pi*angle(VNR1(3,sw1(2))),...
    abs(VNR3(1,sw2(2))),180/pi*angle(VNR3(3,sw2(2))),abs(VNR4(3,sw2(2))),180/pi*angle(VNR4(3,sw2(2))));
fprintf(fid,'\n');

magstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f & %0.4f & %0.4f \\\\\n',...
    (abs(VNR0(1,sw1(1))) - abs(VNR0(1,sw1(2)))),(abs(VNR1(1,sw1(1))) - abs(VNR1(1,sw1(2)))),...
    (abs(VNR3(1,sw2(1))) - abs(VNR3(1,sw2(2)))),(abs(VNR4(1,sw2(1))) - abs(VNR4(1,sw2(2)))));
magstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f & %0.4f & %0.4f \\\\\n',...
    (abs(VNR0(2,sw1(1))) - abs(VNR0(2,sw1(2)))),(abs(VNR1(2,sw1(1))) - abs(VNR1(2,sw1(2)))),...
    (abs(VNR3(2,sw2(1))) - abs(VNR3(2,sw2(2)))),(abs(VNR4(2,sw2(1))) - abs(VNR4(2,sw2(2)))));
magstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f & %0.4f & %0.4f \\\\\n',...
    (abs(VNR0(3,sw1(1))) - abs(VNR0(3,sw1(2)))),(abs(VNR1(3,sw1(1))) - abs(VNR1(3,sw1(2)))),...
    (abs(VNR3(3,sw2(1))) - abs(VNR3(3,sw2(2)))),(abs(VNR4(3,sw2(1))) - abs(VNR4(3,sw2(2)))));
fprintf(fid,'\n');

angstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f & %0.4f & %0.4f \\\\\n',...
    180/pi*(angle(VNR0(1,sw1(1))) - angle(VNR0(1,sw1(2)))),180/pi*(angle(VNR1(1,sw1(1))) - angle(VNR1(1,sw1(2)))),...
    180/pi*(angle(VNR3(1,sw2(1))) - angle(VNR3(1,sw2(2)))),180/pi*(angle(VNR4(1,sw2(1))) - angle(VNR4(1,sw2(2)))));
angstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f & %0.4f & %0.4f \\\\\n',...
    180/pi*(angle(VNR0(2,sw1(1))) - angle(VNR0(2,sw1(2)))),180/pi*(angle(VNR1(2,sw1(1))) - angle(VNR1(2,sw1(2)))),...
    180/pi*(angle(VNR3(2,sw2(1))) - angle(VNR3(2,sw2(2)))),180/pi*(angle(VNR4(2,sw2(1))) - angle(VNR4(2,sw2(2)))));
angstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f & %0.4f & %0.4f \\\\\n',...
    180/pi*(angle(VNR0(3,sw1(1))) - angle(VNR0(3,sw1(2)))),180/pi*(angle(VNR1(3,sw1(1))) - angle(VNR1(3,sw1(2)))),...
    180/pi*(angle(VNR3(3,sw2(1))) - angle(VNR3(3,sw2(2)))),180/pi*(angle(VNR4(3,sw2(1))) - angle(VNR4(3,sw2(2)))));
fprintf(fid,'\n');

powstrA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    SS0(1),SS0(1)/1j,SS1(1),SS1(1)/1j,SS3(1),SS3(1)/1j,SS4(1),SS4(1)/1j);
powstrB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    SS0(2),SS0(2)/1j,SS1(2),SS1(2)/1j,SS3(2),SS3(2)/1j,SS4(2),SS4(2)/1j);
powstrC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',...
    SS0(3),SS0(3)/1j,SS1(3),SS1(3)/1j,SS3(3),SS3(3)/1j,SS4(3),SS4(3)/1j);
fprintf(fid,'\n');

for k1 = 1:network3.nodes.nnode
    if sum(network3.cons.wmaxpu(:,k1)) ~= 0

        fprintf(fid,char(network3.nodes.nodelist(k1)));
        for ph = 1:3
            if network3.cons.wmaxpu(ph,k1) == 0
                fprintf(fid,' & -');                
            else                
                fprintf(fid,' & %0.4f + j%0.4f',real(wopt1(ph,k1)),imag(wopt1(ph,k1)));
            end
        end
        fprintf(fid,' \\\\\n');
    end

end

fclose(fid);

%%

% fid = fopen(['~/Desktop/temp/switch/' name1 '.txt'],'w');
% 
% Vtn1strA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(1,tn1),VNR0(1,tn1)/1j,VNR1(1,tn1),VNR1(1,tn1)/1j,VNR2(1,tn1),VNR2(1,tn1)/1j);
% Vtn1strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(2,tn1),VNR0(2,tn1)/1j,VNR1(2,tn1),VNR1(2,tn1)/1j,VNR2(2,tn1),VNR2(2,tn1)/1j);
% Vtn1strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(3,tn1),VNR0(3,tn1)/1j,VNR1(3,tn1),VNR1(3,tn1)/1j,VNR2(3,tn1),VNR2(3,tn1)/1j);
% fprintf(fid,'\n');
% 
% Vtn1strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
%     abs(VNR0(1,tn1)),180/pi*angle(VNR0(1,tn1)),abs(VNR1(1,tn1)),180/pi*angle(VNR1(1,tn1)),abs(VNR2(1,tn1)),180/pi*angle(VNR2(1,tn1)));
% Vtn1strA = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
%     abs(VNR0(2,tn1)),180/pi*angle(VNR0(2,tn1)),abs(VNR1(2,tn1)),180/pi*angle(VNR1(2,tn1)),abs(VNR2(2,tn1)),180/pi*angle(VNR2(2,tn1)));
% Vtn1strA = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
%     abs(VNR0(3,tn1)),180/pi*angle(VNR0(3,tn1)),abs(VNR1(3,tn1)),180/pi*angle(VNR1(3,tn1)),abs(VNR2(3,tn1)),180/pi*angle(VNR2(3,tn1)));
% fprintf(fid,'\n');
% 
% Vtn2strA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(1,tn2),VNR0(1,tn2)/1j,VNR1(1,tn2),VNR1(1,tn2)/1j,VNR2(1,tn2),VNR2(1,tn2)/1j);
% Vtn2strB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(2,tn2),VNR0(2,tn2)/1j,VNR1(2,tn2),VNR1(2,tn2)/1j,VNR2(2,tn2),VNR2(2,tn2)/1j);
% Vtn2strC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',VNR0(3,tn2),VNR0(3,tn2)/1j,VNR1(3,tn2),VNR1(3,tn2)/1j,VNR2(3,tn2),VNR2(3,tn2)/1j);
% fprintf(fid,'\n');
% 
% Vtn2strA = fprintf(fid,'& $a$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
%     abs(VNR0(1,tn2)),180/pi*angle(VNR0(1,tn2)),abs(VNR1(1,tn2)),180/pi*angle(VNR1(1,tn2)),abs(VNR2(1,tn2)),180/pi*angle(VNR2(1,tn2)));
% Vtn2strA = fprintf(fid,'& $b$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
%     abs(VNR0(2,tn2)),180/pi*angle(VNR0(2,tn2)),abs(VNR1(2,tn2)),180/pi*angle(VNR1(2,tn2)),abs(VNR2(2,tn2)),180/pi*angle(VNR2(2,tn2)));
% Vtn2strA = fprintf(fid,'& $c$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ & $%0.4f \\angle {%0.4f}^{\\circ}$ \\\\\n',...
%     abs(VNR0(3,tn2)),180/pi*angle(VNR0(3,tn2)),abs(VNR1(3,tn2)),180/pi*angle(VNR1(3,tn2)),abs(VNR2(3,tn2)),180/pi*angle(VNR2(3,tn2)));
% fprintf(fid,'\n');
% 
% magstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f & %0.4f \\\\\n',(abs(VNR0(1,tn1)) - abs(VNR0(1,tn2))),(abs(VNR1(1,tn1)) - abs(VNR1(1,tn2))),(abs(VNR2(1,tn1)) - abs(VNR2(1,tn2))));
% magstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f & %0.4f \\\\\n',(abs(VNR0(2,tn1)) - abs(VNR0(2,tn2))),(abs(VNR1(2,tn1)) - abs(VNR1(2,tn2))),(abs(VNR2(2,tn1)) - abs(VNR2(2,tn2))));
% magstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f & %0.4f \\\\\n',(abs(VNR0(3,tn1)) - abs(VNR0(3,tn2))),(abs(VNR1(3,tn1)) - abs(VNR1(3,tn2))),(abs(VNR2(3,tn1)) - abs(VNR2(3,tn2))));
% fprintf(fid,'\n');
% 
% angstrA = fprintf(fid,'& $a$ & %0.4f & %0.4f & %0.4f \\\\\n',180/pi*(angle(VNR0(1,tn1)) - angle(VNR0(1,tn2))),180/pi*(angle(VNR1(1,tn1)) - angle(VNR1(1,tn2))),180/pi*(angle(VNR2(1,tn1)) - angle(VNR2(1,tn2))));
% angstrB = fprintf(fid,'& $b$ & %0.4f & %0.4f & %0.4f \\\\\n',180/pi*(angle(VNR0(2,tn1)) - angle(VNR0(2,tn2))),180/pi*(angle(VNR1(2,tn1)) - angle(VNR1(2,tn2))),180/pi*(angle(VNR2(2,tn1)) - angle(VNR2(2,tn2))));
% angstrC = fprintf(fid,'& $c$ & %0.4f & %0.4f & %0.4f \\\\\n',180/pi*(angle(VNR0(3,tn1)) - angle(VNR0(3,tn2))),180/pi*(angle(VNR1(3,tn1)) - angle(VNR1(3,tn2))),180/pi*(angle(VNR2(3,tn1)) - angle(VNR2(3,tn2))));
% fprintf(fid,'\n');specials concrete jungle
% 
% powstrA = fprintf(fid,'& $a$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',SNR0(1),SNR0(1)/1j,SNR1(1),SNR1(1)/1j,SNR2(1),SNR2(1)/1j);
% powstrB = fprintf(fid,'& $b$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',SNR0(2),SNR0(2)/1j,SNR1(2),SNR1(2)/1j,SNR2(2),SNR2(2)/1j);
% powstrC = fprintf(fid,'& $c$ & %0.4f + j%0.4f & %0.4f + j%0.4f & %0.4f + j%0.4f \\\\\n',SNR0(3),SNR0(3)/1j,SNR1(3),SNR1(3)/1j,SNR2(3),SNR2(3)/1j);
% fprintf(fid,'\n');
% 
% for k1 = 1:nodes1.nnode
%     if sum(controllers1.wmaxpu(:,k1)) > 0
% 
%         fprintf(fid,char(nodes1.nodelist(k1)));
%         for ph = 1:3
%             if controllers1.wmaxpu(ph,k1) == 0
%                 fprintf(fid,' & 0');                
%             else                
%                 fprintf(fid,' & %0.4f + j%0.4f',real(wopt2(ph,k1)),imag(wopt2(ph,k1)));
%             end
%         end
%         fprintf(fid,' \\\\\n');
%     end
% 
% end
% 
% fclose(fid);