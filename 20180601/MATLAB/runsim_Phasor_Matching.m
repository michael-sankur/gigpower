% Michael Sankur - msankur@lbl.gov
% 2018.06.01

clc, clear all, close all

path(path,genpath('/home/michael/Dropbox/Unbalanced LinDistflow/20180601/MATLAB/'));

%% Load feeder

% Change the file path as necessary
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180601/NETWORKS/';
% Change file name as necessary
fn = 'ieee_13node_mesh_open.txt';

[network1] = network_mapper_function(fp, fn);

%% Network paramaters

nnode = network1.nodes.nnode;
nline = network1.lines.nline;

%% Load parameters

network1.loads.aPQ = 0.75*ones(3,nnode).*network1.nodes.PH;
network1.loads.aI = 0.10*ones(3,nnode).*network1.nodes.PH;
network1.loads.aZ = 0.15*ones(3,nnode).*network1.nodes.PH;

network1.loads.spu(:,3:15) = 0.75*network1.loads.spu(:,3:15);
network1.loads.spu(:,16:28) = 1.5*network1.loads.spu(:,16:28);

%% Capacitor parameters

network1.caps.cappu = 1*network1.caps.cappu;

%% Controller parameters

network1.cons.wmaxpu = 0.5*network1.cons.wmaxpu;

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

tn1 = 12;
tn2 = 25;

%% Run CVX

[nvar1, Aineq1, bineq1, Aeq1, beq1] = create_linear_constraints_base(network1,slacknode,Vslack,sim);

cvx_begin quiet
    expressions Zmag1 Zang1 Zw1;
    variable X1(3*nvar1);
    Zmag1 = (X1(tn1) - X1(tn2))^2 + (X1(nvar1+tn1) - X1(nvar1+tn2))^2 + (X1(2*nvar1+tn1) - X1(2*nvar1+tn2))^2;
    Zang1 = (X1(nnode+tn1) - X1(nnode+tn2))^2 + (X1(nvar1+nnode+tn1) - X1(nvar1+nnode+tn2))^2 ...
        + (X1(2*nvar1+nnode+tn1) - X1(2*nvar1+nnode+tn2))^2;
    for ph = 1:3
        for k1 = 1:nnode
            if network1.cons.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar1 + nnode + nnode + nline + nline + k1;
                kv = (ph-1)*nvar1 + nnode + nnode + nline + nline + nnode + k1;
                Zw1 = Zw1 + X1(ku)^2 + X1(kv)^2;
            end
        end
    end
    minimize(1000*Zmag1 + 0*Zang1 + 1*Zw1)
    subject to
    Aeq1 * X1 == beq1
    Aineq1 * X1 <= bineq1
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

[CVXopt1, Eopt1, Topt1, Vopt1, Iopt1, Sopt1, wopt1, dopt1, sopt1] = parse_CVX_output(X1,nvar1,network1);

%% Run CVX 

[nvar2, Aineq2, bineq2, Aeq2, beq2] = create_linear_constraints_base(network1,slacknode,Vslack,sim);

cvx_begin quiet
    expressions Zmag2 Zang2 Zw2;
    variable X2(3*nvar1);
    Zmag2 = (X2(tn1) - X2(tn2))^2 + (X2(nvar2+tn1) - X2(nvar2+tn2))^2 + (X2(2*nvar2+tn1) - X2(2*nvar2+tn2))^2;
    Zang2 = (X2(nnode+tn1) - X2(nnode+tn2))^2 + (X2(nvar2+nnode+tn1) - X2(nvar2+nnode+tn2))^2 ...
        + (X2(2*nvar2+nnode+tn1) - X2(2*nvar2+nnode+tn2))^2;
    for ph = 1:3
        for k1 = 1:nnode
            if network1.cons.wmaxpu(ph,k1) > 0
                ku = (ph-1)*nvar2 + nnode + nnode +  nline + nline + k1;
                kv = (ph-1)*nvar2 + nnode + nnode + nline + nline + nnode + k1;
                Zw2 = Zw2 + X2(ku)^2 + X2(kv)^2;
            end
        end
    end
    minimize(1000*Zmag2 + 1000*Zang2 + 1*Zw2)
    subject to
    Aeq2 * X2 == beq2
    Aineq2 * X2 <= bineq2
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

[CVXopt2, Eopt2, Topt2, Vopt2, Iopt2, Sopt2, wopt2, dopt2, sopt2] = parse_CVX_output(X2,nvar2,network1);

%%

network1.cons.wpu = zeros(3,nnode);
[~,VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, iter0] = NR3(network1,slacknode,Vslack,[],[]);

network1.cons.wpu = wopt1;
[~,VNR1, INR1, STXNR1, SRXNR1, iNR1, sNR1, iter1] = NR3(network1,slacknode,Vslack,Vopt1,Iopt1);

network1.cons.wpu = wopt2;
[~,VNR2, INR2, STXNR2, SRXNR2, iNR2, sNR2, iter2] = NR3(network1,slacknode,Vslack,Vopt2,Iopt2);

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
PD0 = VNR0(:,tn1) - VNR0(:,tn2)
MD0 = (abs(VNR0(:,tn1)) - abs(VNR0(:,tn2)))
AD0 = 180/pi*(angle(VNR0(:,tn1)) - angle(VNR0(:,tn2)))
SF0 = VNR0(:,tn1).*conj(swFYpu*(VNR0(:,tn1) - VNR0(:,tn2)))

disp 'Open, Mag Only Control'
PD1 = VNR1(:,tn1) - VNR1(:,tn2)
MD1 = (abs(VNR1(:,tn1)) - abs(VNR1(:,tn2)))
AD1 = 180/pi*(angle(VNR1(:,tn1)) - angle(VNR1(:,tn2)))
SF1 = VNR1(:,tn1).*conj(swFYpu*(VNR1(:,tn1) - VNR1(:,tn2)))

disp 'Open, Phasor Control'
PD2 = VNR2(:,tn1) - VNR2(:,tn2)
MD2 = (abs(VNR2(:,tn1)) - abs(VNR2(:,tn2)))
AD2 = 180/pi*(angle(VNR2(:,tn1)) - angle(VNR2(:,tn2)))
SF2 = VNR2(:,tn1).*conj(swFYpu*(VNR2(:,tn1) - VNR2(:,tn2)))

disp 'Open, No Control - Open, Mag Only Control - Open, Phasor Control'
disp 'Phasor Difference'
[PD0 PD1 PD2]
disp 'Magnitude Difference'
[MD0 MD1 MD2]
disp 'Angle Difference'
[AD0 AD1 AD2]
disp 'Power Flow'
[SF0 SF1 SF2]
