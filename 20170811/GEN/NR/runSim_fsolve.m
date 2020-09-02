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

% name = '4_node_fullphase';
% name = '4_node_multiphase_tx';
name = 'ieee_13node_new';
% name = '6node_fullphase_new';

fp = ['C:\Users\Michael\Desktop\reconfiguration\Feeders\'];
fn = [name '.txt'];

[feeder, nodes, lines, configs, loads, caps, controllers] = feeder_mapper_TX_function_20170215(name, fp, fn);


%% Feeder paramaters

nnode = nodes.nnode;
nline = lines.nline;

%% Load parameters

loads.aPQ = 1*ones(3,nnode).*nodes.PH;
loads.aI = zeros(3,nnode);
loads.aZ = 0*ones(3,nnode).*nodes.PH;

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



%%

sim.Vfbs = [];
sim.Lfbs = [];
sim.Hfbs = [];
sim.rho = 1;

OPTSOL = Solver_LinDist3Flow_TX_minP0_20160603_magangle(feeder, nodes, lines, configs, loads, caps, controllers, sim)

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

%% Base Case

if radflag == 1    
    FBS0 = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);    
end

createfsolvefunction(feeder,nodes,lines,configs,loads,caps,controllers,'base');

[XFS,fval,ef,output,JFS] = fsolve(@ftemp3phasebase,IG,optimset('MaxIter',1000000,'MaxFunEvals',1000000,'TolFun',1e-6));

for k1 = 2:2:3*2*nnode
    VFS(k1/2) = XFS(k1-1) + 1j*XFS(k1);
end
VFS0 = [VFS(1:nnode);
    VFS(nnode+1:2*nnode);
    VFS(2*nnode+1:end)]

for k1 = 2:2:3*2*nline
    IFS(k1/2) = XFS(3*2*nnode + k1-1) + 1j*XFS(3*2*nnode + k1);
end
IFS0 = [IFS(1:nline);
    IFS(nline+1:2*nline);
    IFS(2*nline+1:end)];

%% Control Case

controllers.wpu = wopt;

if radflag == 1    
    FBS1 = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);    
end

createfsolvefunction(feeder,nodes,lines,configs,loads,caps,controllers,'control');

[XFS,fval,ef,output,JFS] = fsolve(@ftemp3phasecontrol,IG,optimset('MaxIter',1000000,'MaxFunEvals',1000000,'TolFun',1e-6));

for k1 = 2:2:3*2*nnode
    VFS(k1/2) = XFS(k1-1) + 1j*XFS(k1);
end
VFS1 = [VFS(1:nnode);
    VFS(nnode+1:2*nnode);
    VFS(2*nnode+1:end)]

for k1 = 2:2:3*2*nline
    IFS(k1/2) = XFS(3*2*nnode + k1-1) + 1j*XFS(3*2*nnode + k1);
end
IFS1 = [IFS(1:nline);
    IFS(nline+1:2*nline);
    IFS(2*nline+1:end)];
