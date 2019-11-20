% Michael Sankur - msankur@lbl.gov
% 2018.06.01

clc, clear all, close all

path(path,genpath('/home/michael/Dropbox/Unbalanced LinDistflow/20180601/MATLAB/'));

%% Load feeder

% Change the file path as necessary
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180601/NETWORKS/';
% Change file name as necessary
% fn = '05node_fullphase_radial.txt';
% fn = 'ieee_13node_balance.txt';
fn = 'ieee_13node_mesh_open.txt';

[network1] = network_mapper_function(fp, fn);

%% Network paramaters

nnode = network1.nodes.nnode;
nline = network1.lines.nline;

%% Load parameters

network1.loads.aPQ = 0.85*ones(3,nnode).*network1.nodes.PH;
network1.loads.aI = 0.0*ones(3,nnode).*network1.nodes.PH;
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

% [NRRES0, VNR0, INR0, STXNR0, SRXNR0] = NR3(network1,slacknode,Vslack,[],[],1e-9);

%%

Vcvx = ones(3,nnode);
Vpf = zeros(3,nnode);

magerr = 1e99;
angerr = 1e99;


iterflag = 1;
itercount = 0;

while iterflag == 1 % max(max(abs(Vcvx - Vpf))) >= 1e-6
    
    [nvar1, Aineq, bineq, Aeq, beq] = create_linear_constraints_sqp(network1,slacknode,Vslack,sim);

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
        Aineq * X1 <= bineq
        Aeq * X1 == beq
    %     for ph = 1:3
    %         for k1 = 3:nline
    %             kP = (ph-1)*nvar1 + nnode + nnode + k1;
    %             kQ = (ph-1)*nvar1 + nnode + nnode + nline + k1;
    %             norm(X([kP kQ]),2) <= 0.2
    %         end
    %     end
        for ph = 1:3
            for k1 = 1:nnode
                ku = (ph-1)*nvar1 + nnode + nnode + nline + nline + k1;
                kv = (ph-1)*nvar1 + nnode + nnode + nline + nline + nnode + k1;
                norm(X1([ku kv]),2) <= 1.0*network1.cons.wmaxpu(ph,k1)
            end
        end
    cvx_end

    [CVXopt1, Eopt1, Topt1, Vopt1, Iopt1, Sopt1, wopt1, dopt1, sopt1] = parse_CVX_output(X1,nvar1,network1);
    
    Vcvx = Vopt1;
    
    
    network1.cons.wpu = 1*wopt1;

    [NRRES1, VNR1, INR1, STXNR1, SRXNR1] = NR3(network1,slacknode,Vslack,Vopt1,Iopt1,1e-9);
    
    Vpf = VNR1;
    
    sim.VNR = VNR1;
    sim.Lmn = zeros(3,nline);
    sim.Hmn = zeros(3,nline);
    for k2 = 1:nline
        sim.Lmn(:,k2) = (network1.lines.FZpu(:,:,k2)*INR1(:,k2)).*conj(INR1(:,k2));
        sim.Hmn(:,k2) = (network1.lines.FZpu(:,:,k2)*INR1(:,k2)).*conj(network1.lines.FZpu(:,:,k2)*INR1(:,k2));        
    end
    
    Vdif = Vcvx - Vpf;
    
    sim.VNR;
    sim.Lmn;
    sim.Hmn;
    
    disp('~~~~~~~~~~~~~~/\~~~~~~~~~~~~~~~~~~/\~~~~~~~~~~~~~~')
%     disp('~~~~~~~~~~~~~/\/\~~~~~~~~~~~~~~~~/\/\~~~~~~~~~~~~~')
    disp('~~~~~~~~~~~~/\/\/\~~~~~~~~~~~~~~/\/\/\~~~~~~~~~~~~')
%     disp('~~~~~~~~~~~/\/\/\/\~~~~~~~~~~~~/\/\/\/\~~~~~~~~~~~')
    disp('~~~~~~~~~~/\/\/\/\/\~~~~~~~~~~/\/\/\/\/\~~~~~~~~~~')
    disp('~~~~~~~~~~\/\/\/\/\/~~~~~~~~~~\/\/\/\/\/~~~~~~~~~~')
%     disp('~~~~~~~~~~~\/\/\/\/~~~~~~~~~~~~\/\/\/\/~~~~~~~~~~~')
    disp('~~~~~~~~~~~~\/\/\/~~~~~~~~~~~~~~\/\/\/~~~~~~~~~~~~')
%     disp('~~~~~~~~~~~~~\/\/~~~~~~~~~~~~~~~~\/\/~~~~~~~~~~~~~')
    disp('~~~~~~~~~~~~~~\/~~~~~~~~~~~~~~~~~~\/~~~~~~~~~~~~~~')
    
    magerr = max(max(abs(abs(VNR1) - abs(Vopt1))))
%     magerr = max(max(abs(abs(VNR1 - Vopt))))
    
    angerr = max(max(abs(180/pi*(angle(VNR1) - angle(Vopt1)))))
%     angerr = max(max(abs(180/pi*angle(VNR1 - Vopt))))
    
%     Sopt - STXNR1
%     
%     Sopt - SRXNR1

    itercount = itercount + 1;
    
    if magerr <= 1e-6 && angerr <= 1e-3
        disp('Convergence')
        iterflag = 0;
    end
    if itercount >= 10
        disp('Max Iterations')
        iterflag = 0;
    end
    
%     pause    
    
    
end

%%

Vcvx = ones(3,nnode);
Vpf = zeros(3,nnode);

magerr = 1e99;
angerr = 1e99;

iterflag = 1;
itercount = 0;

while iterflag == 1 % max(max(abs(Vcvx - Vpf))) >= 1e-6
    
    [nvar2, Aineq, bineq, Aeq, beq] = create_linear_constraints_sqp(network1,slacknode,Vslack,sim);

    cvx_begin quiet
        expressions Zmag2 Zang2 Zw2;
        variable X2(3*nvar2);
        Zmag2 = (X2(tn1) - X2(tn2))^2 + (X2(nvar2+tn1) - X2(nvar2+tn2))^2 + (X2(2*nvar2+tn1) - X2(2*nvar2+tn2))^2;
        Zang2 = (X2(nnode+tn1) - X2(nnode+tn2))^2 + (X2(nvar2+nnode+tn1) - X2(nvar2+nnode+tn2))^2 ...
            + (X2(2*nvar2+nnode+tn1) - X2(2*nvar2+nnode+tn2))^2;
        for ph = 1:3
            for k1 = 1:nnode
                if network1.cons.wmaxpu(ph,k1) > 0
                    ku = (ph-1)*nvar2 + nnode + nnode + nline + nline + k1;
                    kv = (ph-1)*nvar2 + nnode + nnode + nline + nline + nnode + k1;
                    Zw2 = Zw2 + X2(ku)^2 + X2(kv)^2;
                end
            end
        end
        minimize(1000*Zmag2 + 1000*Zang2 + 1*Zw2)
        subject to
        Aineq * X2 <= bineq
        Aeq * X2 == beq
    %     for ph = 1:3
    %         for k1 = 3:nline
    %             kP = (ph-1)*nvar2 + nnode + nnode + k1;
    %             kQ = (ph-1)*nvar2 + nnode + nnode + nline + k1;
    %             norm(X([kP kQ]),2) <= 0.2
    %         end
    %     end
        for ph = 1:3
            for k1 = 1:nnode
                ku = (ph-1)*nvar2 + nnode + nnode + nline + nline + k1;
                kv = (ph-1)*nvar2 + nnode + nnode + nline + nline + nnode + k1;
                norm(X2([ku kv]),2) <= 1.0*network1.cons.wmaxpu(ph,k1)
            end
        end
    cvx_end

    [CVXopt2, Eopt2, Topt2, Vopt2, Iopt2, Sopt2, wopt2, dopt2, sopt2] = parse_CVX_output(X2,nvar2,network1);
    
    Vcvx = Vopt2;
    
    network1.cons.wpu = 1*wopt2;

    [NRRES2, VNR2, INR2, STXNR2, SRXNR2] = NR3(network1,slacknode,Vslack,Vopt2,Iopt2,1e-9);
    
    Vpf = VNR2;
    
    sim.VNR = VNR2;
    sim.Lmn = zeros(3,nline);
    sim.Hmn = zeros(3,nline);
    for k2 = 1:nline
        sim.Lmn(:,k2) = (network1.lines.FZpu(:,:,k2)*INR2(:,k2)).*conj(INR2(:,k2));
        sim.Hmn(:,k2) = (network1.lines.FZpu(:,:,k2)*INR2(:,k2)).*conj(network1.lines.FZpu(:,:,k2)*INR2(:,k2));        
    end
    
    Vdif = Vcvx - Vpf;
    
    sim.VNR;
    sim.Lmn;
    sim.Hmn;
    
    disp('~~~~~~~~~~~~~~/\~~~~~~~~~~~~~~~~~~/\~~~~~~~~~~~~~~')
%     disp('~~~~~~~~~~~~~/\/\~~~~~~~~~~~~~~~~/\/\~~~~~~~~~~~~~')
    disp('~~~~~~~~~~~~/\/\/\~~~~~~~~~~~~~~/\/\/\~~~~~~~~~~~~')
%     disp('~~~~~~~~~~~/\/\/\/\~~~~~~~~~~~~/\/\/\/\~~~~~~~~~~~')
    disp('~~~~~~~~~~/\/\/\/\/\~~~~~~~~~~/\/\/\/\/\~~~~~~~~~~')
    disp('~~~~~~~~~~\/\/\/\/\/~~~~~~~~~~\/\/\/\/\/~~~~~~~~~~')
%     disp('~~~~~~~~~~~\/\/\/\/~~~~~~~~~~~~\/\/\/\/~~~~~~~~~~~')
    disp('~~~~~~~~~~~~\/\/\/~~~~~~~~~~~~~~\/\/\/~~~~~~~~~~~~')
%     disp('~~~~~~~~~~~~~\/\/~~~~~~~~~~~~~~~~\/\/~~~~~~~~~~~~~')
    disp('~~~~~~~~~~~~~~\/~~~~~~~~~~~~~~~~~~\/~~~~~~~~~~~~~~')
    
    magerr = max(max(abs(abs(VNR2) - abs(Vopt2))))
%     magerr = max(max(abs(abs(VNR1 - Vopt))))
    
    angerr = max(max(abs(180/pi*(angle(VNR2) - angle(Vopt2)))))
%     angerr = max(max(abs(180/pi*angle(VNR1 - Vopt))))
    
%     Sopt - STXNR1
%     
%     Sopt - SRXNR1

    itercount = itercount + 1;
    
    if magerr <= 1e-6 && angerr <= 1e-4
        disp('Convergence')
        iterflag = 0;
    end
    if itercount >= 10
        disp('Max Iterations')
        iterflag = 0;
    end
    
%     pause   
    
end

%%

network1.cons.wpu = zeros(3,nnode);
[~,VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, iter0] = NR3(network1,slacknode,Vslack,[],[],1e-9);

network1.cons.wpu = wopt1;
[~,VNR1, INR1, STXNR1, SRXNR1, iNR1, sNR1, iter1] = NR3(network1,slacknode,Vslack,Vopt1,Iopt1,1e-9);

network1.cons.wpu = wopt2;
[~,VNR2, INR2, STXNR2, SRXNR2, iNR2, sNR2, iter2] = NR3(network1,slacknode,Vslack,Vopt2,Iopt2,1e-9);

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

% disp 'Open, No Control'
PD0 = VNR0(:,tn1) - VNR0(:,tn2);
MD0 = (abs(VNR0(:,tn1)) - abs(VNR0(:,tn2)));
AD0 = 180/pi*(angle(VNR0(:,tn1)) - angle(VNR0(:,tn2)));
SF0 = VNR0(:,tn1).*conj(swFYpu*(VNR0(:,tn1) - VNR0(:,tn2)));

% disp 'Open, Mag Only Control'
PD1 = VNR1(:,tn1) - VNR1(:,tn2);
MD1 = (abs(VNR1(:,tn1)) - abs(VNR1(:,tn2)));
AD1 = 180/pi*(angle(VNR1(:,tn1)) - angle(VNR1(:,tn2)));
SF1 = VNR1(:,tn1).*conj(swFYpu*(VNR1(:,tn1) - VNR1(:,tn2)));

% disp 'Open, Phasor Control'
PD2 = VNR2(:,tn1) - VNR2(:,tn2);
MD2 = (abs(VNR2(:,tn1)) - abs(VNR2(:,tn2)));
AD2 = 180/pi*(angle(VNR2(:,tn1)) - angle(VNR2(:,tn2)));
SF2 = VNR2(:,tn1).*conj(swFYpu*(VNR2(:,tn1) - VNR2(:,tn2)));

disp 'Open, No Control - Open, Mag Only Control - Open, Phasor Control'
disp 'Phasor Difference'
[PD0 PD1 PD2]
disp 'Magnitude Difference'
[MD0 MD1 MD2]
disp 'Angle Difference'
[AD0 AD1 AD2]
disp 'Power Flow'
[SF0 SF1 SF2]

%%

disp(' ')

disp('xxxxxxxxxxxxxxxxxxxxxxx  ʌ  xxxxxxxxxxxxxxxxxxxxxxx')
disp('xxxxxxxxxxxxxxxxxxxxxx  / \  xxxxxxxxxxxxxxxxxxxxxx')
disp('xxxxxxxxxxxxxxxxxxxxx  /   \  xxxxxxxxxxxxxxxxxxxxx')
disp('xxxxxxxxxxxxxxxxxxxxx  | | |  xxxxxxxxxxxxxxxxxxxxx')
disp('xxxxxxxxxxxxxxxxxxxxx  \ \/   xxxxxxxxxxxxxxxxxxxxx')
disp('xxxxxxxxxxxxxxxxxxxxx  /\ \   xxxxxxxxxxxxxxxxxxxxx')
disp('xxxxxxxxxxxxxxxxxxxxx  | | |  xxxxxxxxxxxxxxxxxxxxx')
disp('xxxxxxxxxxxxxxxxxxxxx  \   /  xxxxxxxxxxxxxxxxxxxxx')
disp('xxxxxxxxxxxxxxxxxxxxxx  \ /  xxxxxxxxxxxxxxxxxxxxxx')
disp('xxxxxxxxxxxxxxxxxxxxxxx  v  xxxxxxxxxxxxxxxxxxxxxxx')

disp(' ')

disp('------------x------------')
disp('-----------x-x-----------')
disp('----------x---x----------')
disp('---------x--x--x---------')
disp('---------x--x--x---------')
disp('----------x--xx----------')
disp('----------xx--x----------')
disp('---------x--x--x---------')
disp('---------x--x--x---------')
disp('----------x---x----------')
disp('-----------x-x-----------')
disp('------------x------------')

disp(' ')

disp('------------x------------')
disp('----------x---x----------')
disp('--------x-------x--------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('--------x-----x-x--------')
disp('--------x-x-----x--------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('--------x-------x--------')
disp('----------x---x----------')
disp('------------x------------')

disp(' ')

disp('------------x------------')
disp('-----------x-x-----------')
disp('----------x---x----------')
disp('---------x-----x---------')
disp('--------x-------x--------')
disp('-------x---------x-------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('-------x-----x---x-------')
disp('--------x-----x-x--------')
disp('---------x-----x---------')
disp('--------x-x-----x--------')
disp('-------x---x-----x-------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('------x-----x-----x------')
disp('-------x---------x-------')
disp('--------x-------x--------')
disp('---------x-----x---------')
disp('----------x---x----------')
disp('-----------x-x-----------')
disp('------------x------------')

disp(' ')

disp('------------x-----x---x-----x------------')
disp('-----------x-----x-----x-----x-----------')
disp('----------x-----x-------x-----x----------')
disp('---------x-----x---------x-----x---------')
disp('--------x-----x-----x-----x-----x--------')
disp('-------x-----x-----x-x-----x-----x-------')
disp('------x-----x-----x---x-----x-----x------')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')

disp(' ')

disp('--------------------x--------------------')
disp('-------------------x-x-------------------')
disp('-----------------x-----x-----------------')
disp('---------------x----x----x---------------')
disp('-------------x-----x-x-----x-------------')
disp('-----------x-----x-----x-----x-----------')
disp('---------x-----x----x----x-----x---------')
disp('-------x-----x-----x-x-----x-----x-------')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-------x-----x-----x-----x-----x-x-------')
disp('-------x-x-----x-----x-----x-----x-------')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-----x-----x-----x-----x-----x-----x-----')
disp('-------x-----x-----x-x-----x-----x-------')
disp('---------x-----x----x----x-----x---------')
disp('-----------x-----x-----x-----x-----------')
disp('-------------x-----x-x-----x-------------')
disp('---------------x----x----x---------------')
disp('-----------------x-----x-----------------')
disp('-------------------x-x-------------------')
disp('--------------------x--------------------')

disp('            X            ')
disp('           X X           ')
disp('          X   X          ')
disp('         X  X  X         ')
disp('         X  X  X         ')
disp('          X  XX          ')
disp('          XX  X          ')
disp('         X  X  X         ')
disp('         X  X  X         ')
disp('          X   X          ')
disp('           X X           ')
disp('            X            ')

disp(' ')

disp('            X            ')
disp('          X   X          ')
disp('        X       X        ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('        X     X X        ')
disp('        X X     X        ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('        X       X        ')
disp('          X   X          ')
disp('            X            ')

disp(' ')

disp('            X            ')
disp('           X X           ')
disp('          X   X          ')
disp('         X     X         ')
disp('        X       X        ')
disp('       X         X       ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('       X     X   X       ')
disp('        X     X X        ')
disp('         X     X         ')
disp('        X X     X        ')
disp('       X   X     X       ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('      X     X     X      ')
disp('       X         X       ')
disp('        X       X        ')
disp('         X     X         ')
disp('          X   X          ')
disp('           X X           ')
disp('            X            ')

disp(' ')

disp('            X     X   X     X            ')
disp('           X     X     X     X           ')
disp('          X     X       X     X          ')
disp('         X     X         X     X         ')
disp('        X     X     X     X     X        ')
disp('       X     X     X X     X     X       ')
disp('      X     X     X   X     X     X      ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')

disp(' ')

disp('                    X                    ')
disp('                   X X                   ')
disp('                 X     X                 ')
disp('               X    X    X               ')
disp('             X     X X     X             ')
disp('           X     X     X     X           ')
disp('         X     X    X    X     X         ')
disp('       X     X     X X     X     X       ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('       X     X     X     X     X X       ')
disp('       X X     X     X     X     X       ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('     X     X     X     X     X     X     ')
disp('       X     X     X X     X     X       ')
disp('         X     X    X    X     X         ')
disp('           X     X     X     X           ')
disp('             X     X X     X             ')
disp('               X    X    X               ')
disp('                 X     X                 ')
disp('                   X X                   ')
disp('                    X                    ')