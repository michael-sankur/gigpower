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

networkname = 'ieee_37node_balance';

fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180101/NETWORKS/';
fn = [networkname '.txt'];

[network1] = network_mapper_function_20180101(networkname, fp, fn);

[network2] = network_mapper_function_20180101(networkname, fp, fn);

%% Network paramaters

nnode1 = network1.nodes.nnode;
nline1 = network1.lines.nline;

nnode2 = network2.nodes.nnode;
nline2 = network2.lines.nline;

%% Load parameters

network1.loads.aPQ = 0.85*ones(3,nnode1).*network1.nodes.PH;
network1.loads.aI = zeros(3,nnode1);
network1.loads.aZ = 0.15*ones(3,nnode1).*network1.nodes.PH;

network1.loads.spu = 0.8*network1.loads.spu;

network2.loads.aPQ = 0.85*ones(3,nnode2).*network1.nodes.PH;
network2.loads.aI = zeros(3,nnode2);
network2.loads.aZ = 0.15*ones(3,nnode2).*network1.nodes.PH;

network2.loads.spu = 1.25*network2.loads.spu;

%% Capacitor parameters

network1.caps.cappu = 1*network1.caps.cappu;
% network1.caps.cappu(:,3:15) = 0.5*network1.caps.cappu(:,3:15);
% network1.caps.cappu(:,16:28) = 0.5*network1.caps.cappu(:,16:28);

network2.caps.cappu = 1*network2.caps.cappu;

%% Controller parameters

network1.cons.wmaxpu = 0.5*network1.cons.wmaxpu;

network2.cons.wmaxpu = 0.5*network2.cons.wmaxpu;


%% Simulation parameters

Vnom1 = 1*[1;
    1*exp(1j*-120*pi/180);
    1*exp(1j*120*pi/180)];
sim1.Vnom = Vnom1;

sim1.VNR = [];
sim1.Lmn = [];
sim1.Hmn = [];

Vnom2 = 1*[1;
    1*exp(1j*-120*pi/180);
    1*exp(1j*120*pi/180)]*1.01*exp(1j*0.5*pi/180);
sim2.Vnom = Vnom2;

sim2.VNR = [];
sim2.Lmn = [];
sim2.Hmn = [];


swtn = [38 38;
    23 23];

%% Run CVX

[nvar1, Aineq1, bineq1, Aeq1, beq1] = create_CVX_Matrix_CloseSwitch_20180101(network1,sim1);

[nvar2, Aineq2, bineq2, Aeq2, beq2] = create_CVX_Matrix_CloseSwitch_20180101(network2,sim2);

cvx_begin quiet
    expressions Zmag1 Zang1 Zw1 Zw2;
    variable X1(3*nvar1);
    variable X2(3*nvar2);
%     Zmag1 = (X1(tn1) - X1(tn2))^2 + (X1(nvar1+tn1) - X1(nvar1+tn2))^2 + (X1(2*nvar1+tn1) - X1(2*nvar1+tn2))^2;
%     Zang1 = (X1(nnode+tn1) - X1(nnode+tn2))^2 + (X1(nvar1+nnode+tn1) - X1(nvar1+nnode+tn2))^2 ...
%         + (X1(2*nvar1+nnode+tn1) - X1(2*nvar1+nnode+tn2))^2;
    Zmag1 = 0;
    Zang1 = 0;
    for k1 = 1:size(swtn,1)
        tn1 = swtn(k1,1);
        tn2 = swtn(k1,2);
        Zmag1 = Zmag1 + (X1(tn1) - X2(tn2))^2 + (X1(nvar1+tn1) - X2(nvar2+tn2))^2 + (X1(2*nvar1+tn1) - X2(2*nvar2+tn2))^2;
        Zang1 = Zang1 + (X1(nnode1+tn1) - X2(nnode2+tn2))^2 + (X1(nvar1+nnode1+tn1) - X2(nvar2+nnode2+tn2))^2 ...
            + (X1(2*nvar1+nnode1+tn1) - X2(2*nvar2+nnode2+tn2))^2;        
    end
    for ph = 1:3
        for k1 = 1:nnode1
            if network1.cons.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + k1;
                kv = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + nnode1 + k1;
                Zw1 = Zw1 + X1(ku)^2 + X1(kv)^2;
            end
        end
    end
    for ph = 1:3
        for k1 = 1:nnode2
            if network1.cons.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar2 + nnode2 + nnode2 + nline2 + nline2 + k1;
                kv = (ph-1)*nvar2 + nnode2 + nnode2 + nline2 + nline2 + nnode2 + k1;
                Zw2 = Zw2 + X2(ku)^2 + X2(kv)^2;
            end
        end
    end
    minimize(1000*Zmag1 + 0*Zang1 + 1*Zw2)
    subject to;
    Aeq1 * X1 == beq1;
    Aineq1 * X1 <= bineq1;
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
%             ku = (ph-1)*nvar1 + nnode + nnode +  nline + nline + k1;
%             kv = (ph-1)*nvar1 + nnode + nnode + nline + nline + nnode + k1;
%             norm(X1([ku kv]),2) <= 1.0*cons.wmaxpu(ph,k1)
%         end
%     end
cvx_end

[Eopt11, Dopt11, Vopt11, Iopt11, Sopt11, wopt11, demopt11, sopt11] = parse_CVX_output(X1,nvar1,network1);

[Eopt12, Dopt12, Vopt12, Iopt12, Sopt12, wopt12, demopt12, sopt12] = parse_CVX_output(X2,nvar2,network2);

%% Run CVX 

[nvar1, Aineq1, bineq1, Aeq1, beq1] = create_CVX_Matrix_CloseSwitch_20180101(network1,sim1);

[nvar2, Aineq2, bineq2, Aeq2, beq2] = create_CVX_Matrix_CloseSwitch_20180101(network2,sim2);

cvx_begin quiet
    expressions Zmag1 Zang1 Zw1 Zw2;
    variable X1(3*nvar1);
    variable X2(3*nvar2);
%     Zmag1 = (X1(tn1) - X1(tn2))^2 + (X1(nvar1+tn1) - X1(nvar1+tn2))^2 + (X1(2*nvar1+tn1) - X1(2*nvar1+tn2))^2;
%     Zang1 = (X1(nnode+tn1) - X1(nnode+tn2))^2 + (X1(nvar1+nnode+tn1) - X1(nvar1+nnode+tn2))^2 ...
%         + (X1(2*nvar1+nnode+tn1) - X1(2*nvar1+nnode+tn2))^2;
    Zmag1 = 0;
    Zang1 = 0;
    for k1 = 1:size(swtn,1)
        tn1 = swtn(k1,1);
        tn2 = swtn(k1,2);
        Zmag1 = Zmag1 + (X1(tn1) - X2(tn2))^2 + (X1(nvar1+tn1) - X2(nvar2+tn2))^2 + (X1(2*nvar1+tn1) - X2(2*nvar2+tn2))^2;
        Zang1 = Zang1 + (X1(nnode1+tn1) - X2(nnode2+tn2))^2 + (X1(nvar1+nnode1+tn1) - X2(nvar2+nnode2+tn2))^2 ...
            + (X1(2*nvar1+nnode1+tn1) - X2(2*nvar2+nnode2+tn2))^2;        
    end
    for ph = 1:3
        for k1 = 1:nnode1
            if network1.cons.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + k1;
                kv = (ph-1)*nvar1 + nnode1 + nnode1 + nline1 + nline1 + nnode1 + k1;
                Zw1 = Zw1 + X1(ku)^2 + X1(kv)^2;
            end
        end
    end
    for ph = 1:3
        for k1 = 1:nnode2
            if network1.cons.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar2 + nnode2 + nnode2 + nline2 + nline2 + k1;
                kv = (ph-1)*nvar2 + nnode2 + nnode2 + nline2 + nline2 + nnode2 + k1;
                Zw2 = Zw2 + X2(ku)^2 + X2(kv)^2;
            end
        end
    end
    minimize(1000*Zmag1 + 1000*Zang1 + 1*Zw2)
    subject to;
    Aeq1 * X1 == beq1;
    Aineq1 * X1 <= bineq1;
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
%             ku = (ph-1)*nvar1 + nnode + nnode +  nline + nline + k1;
%             kv = (ph-1)*nvar1 + nnode + nnode + nline + nline + nnode + k1;
%             norm(X1([ku kv]),2) <= 1.0*cons.wmaxpu(ph,k1)
%         end
%     end
cvx_end

[Eopt21, Dopt21, Vopt21, Iopt21, Sopt21, wopt21, demopt21, sopt21] = parse_CVX_output(X1,nvar1,network1);

[Eopt22, Dopt22, Vopt22, Iopt22, Sopt22, wopt22, demopt22, sopt22] = parse_CVX_output(X2,nvar2,network2);

%%

network1.cons.wpu = zeros(3,nnode1);
[VNR01, INR01, STXNR01, SRXNR01, iNR01, sNR01, iter01] = NR3(network1,[],[],1,Vnom1);

network2.cons.wpu = zeros(3,nnode2);
[VNR02, INR02, STXNR02, SRXNR02, iNR02, sNR02, iter02] = NR3(network2,[],[],1,Vnom2);

%%
network1.cons.wpu = wopt11;
[VNR11, INR11, STXNR11, SRXNR11, iNR11, sNR11, iter11] = NR3(network1,Vopt11,Iopt11,1,Vnom1);

network2.cons.wpu = wopt12;
[VNR12, INR12, STXNR12, SRXNR12, iNR12, sNR12, iter12] = NR3(network2,Vopt12,Iopt12,1,Vnom2);

%%
network1.cons.wpu = wopt21;
[VNR21, INR21, STXNR21, SRXNR21, iNR21, sNR21, iter21] = NR3(network1,Vopt21,Iopt21,1,Vnom1);

network2.cons.wpu = wopt22;
[VNR22, INR22, STXNR22, SRXNR22, iNR22, sNR22, iter22] = NR3(network2,Vopt22,Iopt22,1,Vnom2);

%%

% swFYpu = [15.8536 -1j*45.6950, -6.7264 + 1j*16.8943, -3.6834 + 1j*12.6291;
%   -6.7264 + 1j*16.8943, 13.8820 - 1j*43.3005, -1.7486 + 1j*9.6454;
%   -3.6834 + 1j*12.6291, -1.7486 + 1j*9.6454, 12.2761 - 1j*40.8493];
swFYpu = ...
    [15.8536, -6.7264, -3.6834;
    -6.7264, 13.8820, -1.7486;
    -3.6834, -1.7486, 12.2761] + ...
    1j*...
    [-45.6950, 16.8943, 12.6291;
   16.8943, -43.3005, 9.6454;
   12.6291, 9.6454, -40.8493];


disp 'Open, No Control'
for k1 = 1:size(swtn,1)
    tn1 = swtn(k1,1);
    tn2 = swtn(k1,2);
    VD0(:,k1) = VNR01(:,tn1) - VNR02(:,tn2);
    MD0(:,k1) = abs(VNR01(:,tn1)) - abs(VNR02(:,tn2));
    AD0(:,k1) = 180/pi*(angle(VNR01(:,tn1)) - angle(VNR02(:,tn2)));
    SS0(:,k1) = VNR01(:,tn1).*conj(swFYpu*(VNR01(:,tn1) - VNR02(:,tn2)));
end

disp 'Open, Mag Only Control'
for k1 = 1:size(swtn,1)
    tn1 = swtn(k1,1);
    tn2 = swtn(k1,2);
    VD1(:,k1) = VNR11(:,tn1) - VNR12(:,tn2);
    MD1(:,k1) = abs(VNR11(:,tn1)) - abs(VNR12(:,tn2));
    AD1(:,k1) = 180/pi*(angle(VNR11(:,tn1)) - angle(VNR12(:,tn2)));
    SS1(:,k1) = VNR11(:,tn1).*conj(swFYpu*(VNR11(:,tn1) - VNR12(:,tn2)));
end

disp 'Open, Phasor Control'
for k1 = 1:size(swtn,1)
    tn1 = swtn(k1,1);
    tn2 = swtn(k1,2);
    VD2(:,k1) = VNR21(:,tn1) - VNR22(:,tn2);
    MD2(:,k1) = abs(VNR21(:,tn1)) - abs(VNR22(:,tn2));
    AD2(:,k1) = 180/pi*(angle(VNR21(:,tn1)) - angle(VNR22(:,tn2)));
    SS2(:,k1) = VNR21(:,tn1).*conj(swFYpu*(VNR21(:,tn1) - VNR22(:,tn2)));
end