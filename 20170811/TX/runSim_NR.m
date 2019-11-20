% Michael Sankur - msankur@lbl.gov
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

if exist('FBS','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
    path(path,'C:\Users\Michael\Desktop\reconfiguration\20160603');
end

if exist('mxlpsolve','file') == 0
    path(path,'C:\Users\Michael\Desktop\mxlp');
end

if exist('cvx_begin','file') == 0
    cd C:\Users\Michael\Desktop\cvx-w64\cvx
    cvx_setup
end

%% Load feeder

name = '4_node_multiphase_tx';
name = '4_node_fullphase';
% name = '3_node_fullphase';
name = 'ieee_13node_new';
% name = '6node_fullphase_new';

fp = ['C:\Users\Michael\Desktop\reconfiguration\Feeders\'];
fn = [name '.txt'];

[feeder, nodes, lines, configs, loads, caps, controllers] = feeder_mapper_TX_function_20170215(name, fp, fn);


%% Feeder paramaters

nnode = nodes.nnode;
nline = lines.nline;

%% Load parameters

loads.aPQ = 0.85*ones(3,nnode).*nodes.PH;
loads.aI = zeros(3,nnode);
loads.aZ = 0.15*ones(3,nnode).*nodes.PH;

loads.spu = 1.0*loads.spu;

%% Capacitor parameters

caps.cappu = 1*caps.cappu;

%% Controller parameters

controllers.wmaxpu = 1*controllers.wmaxpu;

controllers.wpu = zeros(3,nnode);

%% FBS New

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];

%%

radflag = 1;
for k1 = 1:size(nodes.FM,1)
   
    if sum(nodes.FM(k1,:) == -1) > 1
        radflag = 0;
    end
    
end
radflag


%%

sim.Vfbs = [];
sim.Lfbs = [];
sim.Hfbs = [];
sim.rho = 1;

OPTSOL = Solver_LinDist3Flow_TX_minP0_20160603_magangle(feeder, nodes, lines, configs, loads, caps, controllers, sim)
% OPTSOL = Solver_LinDist3Flow_TX_min36_20160603_magangle(feeder, nodes, lines, configs, loads, caps, controllers, sim)


nvar = OPTSOL.nvar;
Xa = OPTSOL.X(1:nvar);
Xb = OPTSOL.X(nvar+1:2*nvar);
Xc = OPTSOL.X(2*nvar+1:3*nvar);

XX = [Xa Xb Xc]';

Yopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

Yopt
Vopt = sqrt(Yopt)

Dopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

Dopt

Popt = XX(:,1:nline);
XX(:,1:nline) = [];

Qopt = XX(:,1:nline);
XX(:,1:nline) = [];

Sopt = Popt + 1j*Qopt

uopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

vopt = XX(:,1:nnode);
XX(:,1:nnode) = [];
clear XX

wopt = uopt + 1j*vopt

dem = loads.spu.*(loads.aPQ + loads.aZ.*Yopt).*nodes.PH

for k1 = 2:nnode
    sopt(:,k1) = sum(Sopt(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
        sum(Sopt(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2);
end
sopt

%%

Vopt = sqrt(Yopt).*exp(j*pi/180*Dopt)

% FBS = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% Vopt = FBS.V;

for k1 = 1:nline
    
    Iopt(:,k1) = lines.FYpu(:,:,k1)*(Vopt(:,lines.TXnum(k1)) - Vopt(:,lines.RXnum(k1)));
        
end

IG = [];
for ph = 1:3    
    for k1 = 1:nnode    
        IG((ph-1)*2*nnode + 2*k1 - 1) = real(Vopt(ph,k1));
        IG((ph-1)*2*nnode + 2*k1) = imag(Vopt(ph,k1));    
    end
end
for ph = 1:3
    for k1= 1:nline
        IG(3*2*nnode + (ph-1)*2*nline + 2*k1 - 1) = real(Iopt(ph,k1));
        IG(3*2*nnode + (ph-1)*2*nline + 2*k1) = imag(Iopt(ph,k1));
    end
end

%%

XNR0 = IG.';

% XNR0 = [];
% for k1 = 1:nodes.nnode
%     XNR0 = [XNR0;
%         real(Vnom(1));
%         imag(Vnom(1));
%         real(Vnom(2));
%         imag(Vnom(2));
%         real(Vnom(3));
%         imag(Vnom(3))];
% end
% XNR0 = [XNR0; 0.1*ones(6*lines.nline,1)];

FT = 1e99;
iter0 = 0;
while max(abs(FT)) > 1e-9

    FT = computeFTreal(XNR0,feeder,nodes,lines,configs,loads,caps,controllers,2,Vnom);

    JT = computeJTreal(XNR0,feeder,nodes,lines,configs,loads,caps,controllers,2);
    
    if size(JT,1) >= size(JT,2)        
        XNR0 = XNR0 - inv(JT.'*JT)*JT.'*FT        
    end
    
    pause
    
    iter0 = iter0 + 1;
        
end

for k1 = 2:2:3*2*nnode
    VNR0(k1/2) = XNR0(k1-1) + 1j*XNR0(k1);
end
VNR0 = [VNR0(1:nnode);
    VNR0(nnode+1:2*nnode);
    VNR0(2*nnode+1:end)]

for k1 = 2:2:3*2*nline
    INR0(k1/2) = XNR0(3*2*nnode + k1-1) + 1j*XNR0(3*2*nnode + k1);
end
INR0 = [INR0(1:nline);
    INR0(nline+1:2*nline);
    INR0(2*nline+1:end)];

createfsolvefunction(feeder,nodes,lines,configs,loads,caps,controllers,'base');

[XFS0,fval,ef,output,JFS] = fsolve(@ftemp3phasebase,IG,optimset('MaxIter',1000000,'MaxFunEvals',1000000,'TolFun',1e-6));

for k1 = 2:2:3*2*nnode
    VFS0(k1/2) = XFS0(k1-1) + 1j*XFS0(k1);
end
VFS0 = [VFS0(1:nnode);
    VFS0(nnode+1:2*nnode);
    VFS0(2*nnode+1:end)]

for k1 = 2:2:3*2*nline
    IFS0(k1/2) = XFS0(3*2*nnode + k1-1) + 1j*XFS0(3*2*nnode + k1);
end
IFS0 = [IFS0(1:nline);
    IFS0(nline+1:2*nline);
    IFS0(2*nline+1:end)];

if radflag == 1    
    FBS0 = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);    
end

%%

controllers.wpu = wopt;

XNR1 = IG.';
FT = 1e99;
iter1 = 0;
while max(abs(FT)) > 1e-9

    FT = computeFTreal(XNR1,feeder,nodes,lines,configs,loads,caps,controllers,2,Vnom);

    JT = computeJTreal(XNR1,feeder,nodes,lines,configs,loads,caps,controllers,2);
    
    if size(JT,1) >= size(JT,2)        
        XNR1 = XNR1 - inv(JT.'*JT)*JT.'*FT;        
    end
    
    iter1 = iter1 + 1;
    
end

for k1 = 2:2:3*2*nnode
    VNR1(k1/2) = XNR1(k1-1) + 1j*XNR1(k1);
end
VNR1 = [VNR1(1:nnode);
    VNR1(nnode+1:2*nnode);
    VNR1(2*nnode+1:end)]

for k1 = 2:2:3*2*nline
    INR1(k1/2) = XNR1(3*2*nnode + k1-1) + 1j*XNR1(3*2*nnode + k1);
end
INR1 = [INR1(1:nline);
    INR1(nline+1:2*nline);
    INR1(2*nline+1:end)];

createfsolvefunction(feeder,nodes,lines,configs,loads,caps,controllers,'control');

[XFS1,fval,ef,output,JFS] = fsolve(@ftemp3phasecontrol,IG,optimset('MaxIter',1000000,'MaxFunEvals',1000000,'TolFun',1e-6));

for k1 = 2:2:3*2*nnode
    VFS1(k1/2) = XFS1(k1-1) + 1j*XFS1(k1);
end
VFS1 = [VFS1(1:nnode);
    VFS1(nnode+1:2*nnode);
    VFS1(2*nnode+1:end)]

for k1 = 2:2:3*2*nline
    IFS1(k1/2) = XFS1(3*2*nnode + k1-1) + 1j*XFS1(3*2*nnode + k1);
end
IFS1 = [IFS1(1:nline);
    IFS1(nline+1:2*nline);
    IFS1(2*nline+1:end)];

if radflag == 1    
    FBS1 = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);    
end

%%

for k1 = 1:nline
    
    STXNR0(:,k1) = VNR0(:,lines.TXnum(k1)).*conj(INR0(:,k1));
    SRXNR0(:,k1) = VNR0(:,lines.RXnum(k1)).*conj(INR0(:,k1));
    
end

for k1 = 1:nline
    
    STXNR1(:,k1) = VNR1(:,lines.TXnum(k1)).*conj(INR1(:,k1));
    SRXNR1(:,k1) = VNR1(:,lines.RXnum(k1)).*conj(INR1(:,k1));
    
end


