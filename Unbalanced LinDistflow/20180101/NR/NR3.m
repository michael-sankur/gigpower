function [VNR, INR, STXNR, SRXNR, iNR, sNR, iter] = NR3(network,Vopt,Iopt,slackidx,slackVnom,tol)

% Michael Sankur - msankur@lbl.gov
% 2018.01.01

% This function runs a Newton-Raphson algorithm to solve power flow 

% INPUT(S)
% network - struct containing all pertinent the network information,
% including all other structs
% base - struct containing base values
% nodes - struct containing node parameters
% lines - struct containing line parameters
% loads - struct containing load parameters
% caps - struct containing capacitor parameters
% cons - struct containing controller parameters
% vvc - struct containing vvc paramaters
% slackidx - index of slack node
% slackVnom - voltage reference for slack node

% OUTPUT(S)
% VNR - Node voltage in 3 x nnode matrix
% INR - Line current in 3 x nline matrix
% STXNR - Line power at sending (TX) end in 3 x nline matrix
% SRXNR - Line power at receiving (RX) end in 3 x nline matrix
% iNR - Total node current in 3 x nnode matrix - sum of currents from
% loads, capacitors, and controllers
% sNR - Total node power in 3 x nnode matrix - sum of the complex power
% from loads, capacitros and controllers

% Vopt and Iopt are the optimal voltage and current from CVX

% slackidx is the node index of the slack bus, which is assigned a fixed
% voltage reference of slackVnom.

% Voltage and current are separated into their real and imaginary parts
% V_n^phi = A_n^phi + j B_n^phi
% I_n^phi = C_n^phi + j D_n^phi

% Voltage and current vectors for a single phase
% V^phi = [A_1^phi, B_1^phi, A_2^phi, B_2^phi, ... , A_n^phi, B_n^phi]
% I^phi = [C_1^phi, D_1^phi, C_2^phi, D_2^phi, ... , C_n^phi, D_n^phi]

% The NR algorithm variable
% X = [V^a V^b V^c I^a I^b I^c]

if ~exist('tol','var')
  tol=1e-9;
end

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
inmat = nodes.inmat;
outmat = nodes.outmat;

% line paramters
nline = lines.nline;
LPH = lines.PH;
TXnum = lines.TXnum;
RXnum = lines.RXnum;
FZpu = lines.FZpu;
FRpu = lines.FRpu;
FXpu = lines.FXpu;

% load parameters
spu = loads.spu;
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
if isempty(Vopt) == 0
    for ph = 1:3
        for k1 = 1:nnode
            XNR(2*(ph-1)*nnode + 2*k1-1) = real(Vopt(ph,k1));
            XNR(2*(ph-1)*nnode + 2*k1) = imag(Vopt(ph,k1));
        end
    end
% If no initial V is given, use nominal 1 pu balanced voltage profile
elseif isempty(Vopt) == 1
    for ph = 1:3
        for k1 = 1:nnode
            XNR(2*(ph-1)*nnode + 2*k1-1) = real(slackVnom(ph));
            XNR(2*(ph-1)*nnode + 2*k1) = imag(slackVnom(ph));
        end
    end
end

% Initialize current portion of XNR
% If initial I is given (usually from CVX)
if isempty(Iopt) == 0
    for ph = 1:3
        for k1 = 1:nline
            XNR(2*3*nnode + 2*(ph-1)*nline + 2*k1-1) = real(Iopt(ph,k1));
            XNR(2*3*nnode + 2*(ph-1)*nline + 2*k1) = imag(Iopt(ph,k1));
%             2*3*nnode + 2*(ph-1)*nnode + 2*k1
        end
    end
% If no initial I is given, use flat current profile
elseif isempty(Iopt) == 1
    XNR(2*3*nnode+1:end) =  0.01*ones(2*3*nline,1);
end


FT = 1e99;
iter = 0;
while max(abs(FT)) >= tol

    FT = computeFTreal(XNR,network,slackidx,slackVnom);

    JT = computeJTreal(XNR,network,slackidx);
    
%     size(XNR);
%     size(FT);
%     size(JT);
    
    if size(JT,1) >= size(JT,2)
%         size((JT.'*JT)\(JT.'*FT))
        XNR = XNR - (JT.'*JT)\(JT.'*FT);
        eig(JT.'*JT);
    end
       
    iter = iter + 1;
        
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

% STXNR_n^phi = V_m^phi (I_mn^phi)^*
% SRXNR_n^phi = V_n^phi (I_mn^phi)^*
for k1 = 1:nline
    STXNR(:,k1) = VNR(:,TXnum(k1)).*conj(INR(:,k1));
    SRXNR(:,k1) = VNR(:,RXnum(k1)).*conj(INR(:,k1));    
end
STXNR(lines.PH == 0) = 0;
SRXNR(lines.PH == 0) = 0;


% Total node loads
sNR = spu.*(aPQ + aZ.*abs(VNR).^2) - 1j*cappu + wpu + 1j*vvcpu;
sNR(nodes.PH == 0) = 0;
% Total node current
iNR = conj(sNR./VNR);
iNR(nodes.PH == 0) = 0;

end