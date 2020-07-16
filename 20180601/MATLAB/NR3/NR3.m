function [NRRES, VNR, INR, STXNR, SRXNR, iNR, sNR, itercount] = NR3(network, slacknode, Vslack, V0, I0, tol, maxiter)

% Michael Sankur - msankur@lbl.gov
% 2018.06.01

% This function initializes and runs a Newton-Raphson algorithm to solve
% power flow on an unbalanced network. All networks are assumed to be three
% phase, and all output matrices will have three rows representing three
% phases, with non-existing phases padded with zeros

% INPUT(S)
% network - struct containing all pertinent the network information,
% including all other structs
% network.base - struct containing base values
% network.nodes - struct containing node parameters
% network.lines - struct containing line parameters
% network.loads - struct containing load parameters
% network.caps - struct containing capacitor parameters
% network.cons - struct containing controller parameters
% network.vvc - struct containing vvc paramaters
% slacknode - index of slack node
% Vslack - voltage reference for slack node
% V0 - matrix of node voltages for initializing NR algorthm
% I0 - matrix of line currents for initializing NR algorthm
% tol - NR algorithm tolerance
% maxiter - maximum number of iterations for NR algorthm

% OUTPUT(S)
% NRRES - struct containing network states from NR algorthm
% VNR - Node voltage in 3 x nnode matrix
% INR - Line current in 3 x nline matrix
% STXNR - Line power at sending (TX) end in 3 x nline matrix
% SRXNR - Line power at receiving (RX) end in 3 x nline matrix
% iNR - Total node current in 3 x nnode matrix - sum of currents from
% loads, capacitors, and controllers
% sNR - Total node power in 3 x nnode matrix - sum of the complex power
% from loads, capacitros and controllers

% Voltage and current are separated into their real and imaginary parts
% V_n^phi = A_n^phi + j B_n^phi
% I_n^phi = C_n^phi + j D_n^phi

% Voltage and current vectors for a phase phi
% V^phi = [A_1^phi, B_1^phi, A_2^phi, B_2^phi, ... , A_n^phi, B_n^phi]
% I^phi = [C_1^phi, D_1^phi, C_2^phi, D_2^phi, ... , C_n^phi, D_n^phi]

% The NR algorithm variable
% XNR = [V^a V^b V^c I^a I^b I^c]

% network parameters
base = network.base;
nodes = network.nodes;
lines = network.lines;
configs = network.configs;
loads = network.loads;
caps = network.caps;
cons = network.cons;
vvc = network.vvc;

% node parameters
nnode = nodes.nnode;
NPH = nodes.PH;
inlines = nodes.inlines;
innodes = nodes.innodes;
outlines = nodes.outlines;
outnodes = nodes.outnodes;

% line paramters
nline = lines.nline;
LPH = lines.PH;
TXnum = lines.TXnum;
RXnum = lines.RXnum;
FZpu = lines.FZpu;
FRpu = lines.FRpu;
FXpu = lines.FXpu;

% load parameters
spu = loads.spu_nom;
aPQ = loads.aPQ;
aI = loads.aI;
aZ = loads.aZ;

% capacitor paramters
cappu = caps.cappu;

% controller parameters
wpu = cons.wpu;

% vvc parameters
vvcpu = vvc.vvcpu;

% Intialize voltage portion of XNR
XNR = zeros(2*3*(nnode+nline),1);
% If initial V is given (usually from CVX)
if isempty(V0) == 0
    for ph = 1:3
        for k1 = 1:nnode
            XNR(2*(ph-1)*nnode + 2*k1-1) = real(V0(ph,k1));
            XNR(2*(ph-1)*nnode + 2*k1) = imag(V0(ph,k1));
        end
    end
% If no initial V is given, use nominal 1 pu balanced voltage profile
elseif isempty(V0) == 1 || ~exist('V0','var')
    for ph = 1:3
        for k1 = 1:nnode
            XNR(2*(ph-1)*nnode + 2*k1-1) = real(Vslack(ph));
            XNR(2*(ph-1)*nnode + 2*k1) = imag(Vslack(ph));
        end
    end
end

% Initialize current portion of XNR
% If initial I is given (usually from CVX)
if isempty(I0) == 0
    for ph = 1:3
        for k1 = 1:nline
            XNR(2*3*nnode + 2*(ph-1)*nline + 2*k1-1) = real(I0(ph,k1));
            XNR(2*3*nnode + 2*(ph-1)*nline + 2*k1) = imag(I0(ph,k1));
%             2*3*nnode + 2*(ph-1)*nnode + 2*k1
        end
    end
% If no initial I is given, use flat current profile
elseif isempty(I0) == 1 || ~exist('I0','var')
    XNR(2*3*nnode+1:end) =  0.01*ones(2*3*nline,1);
end

% If no tolerance given, set to 1e-9
if ~exist('tol','var')
  tol=1e-9;
end

% If no maximum iterations given, set to 100
if ~exist('maxiter','var')
  maxiter = 100;
end

% Newton-Raphson for loop
FT = 1e99;
itercount = 0;

% Loop until all residuals within tolerance bound
while max(abs(FT)) >= tol && itercount < maxiter

    % Compute residuals
    FT = computeFTreal(XNR,network,slacknode,Vslack);
    
    % Compute derivaties of residuals
    JT = computeJTreal(XNR,network,slacknode);
    
        
    % If Jacobian matrix is tall
    if size(JT,1) >= size(JT,2)
        XNR = XNR - (JT.'*JT)\(JT.'*FT);
        % eig(JT.'*JT);  
    end   
    itercount = itercount + 1;  
end

% remap XNR to VNR, INR, STXNR, SRXNR, iNR, sNR
% VNR = XNR(1:2:2*3*nnode-1).' + 1j*XNR(2:2:2*3*nnode).';
for k1 = 2:2:3*2*nnode
    VNR(k1/2) = XNR(k1-1) + 1j*XNR(k1);
end
VNR = [VNR(1:nnode);
    VNR(nnode+1:2*nnode);
    VNR(2*nnode+1:end)];
VNR(nodes.PH == 0) = 0;
% INR = XNR(2*3*nnode+1:2:2*3*nnode+2*3*nline-1) + 1j*XNR(2*3*nnode+2:2:2*3*nnode+2*3*nline)
for k1 = 2:2:3*2*nline
    INR(k1/2) = XNR(3*2*nnode + k1-1) + 1j*XNR(3*2*nnode + k1);
end
INR = [INR(1:nline);
    INR(nline+1:2*nline);
    INR(2*nline+1:end)];
INR(lines.PH == 0) = 0;

% STXNR_m^phi = V_m^phi (I_mn^phi)^*
% SRXNR_n^phi = V_n^phi (I_mn^phi)^*
for k1 = 1:nline
    STXNR(:,k1) = VNR(:,TXnum(k1)).*conj(INR(:,k1));
    SRXNR(:,k1) = VNR(:,RXnum(k1)).*conj(INR(:,k1));    
end
STXNR(lines.PH == 0) = 0;
SRXNR(lines.PH == 0) = 0;

% Total node loads
sNR = spu.*(aPQ + aI.*abs(VNR) + aZ.*abs(VNR).^2) - 1j*cappu + wpu + 1j*vvcpu;
sNR(nodes.PH == 0) = 0;
% Total node current
iNR = conj(sNR./VNR);
iNR(nodes.PH == 0) = 0;

NRRES.VNR = VNR;
NRRES.INR = INR;
NRRES.STXNR = STXNR;
NRRES.SRXNR = SRXNR;
NRRES.iNR = iNR;
NRRES.sNR = sNR;
NRRES.itercount = itercount;

end