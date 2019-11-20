function [Eopt, Dopt, Vopt, Iopt, Sopt, wopt, demopt, sopt] = parse_CVX_output(X,nvar,network)

% Michael Sankur - msankur@lbl.gov
% 2018.01.01

% This function parses the output of CVX, which is the variable X

% INPUT(S)
% X - optimal solution from CVX
% nvar - number of variables per phase
% network - struct containing all pertinent network information,
% including all other structs
% base - struct containing base values
% nodes - struct containing node parameters
% lines - struct containing line parameters
% loads - struct containing load parameters
% caps - struct containing capacitor parameters
% cons - struct containing controller parameters

% OUTPUT(S)
% Eopt - node squared voltage magnitude in 3 x nnode matrix
% Dopt - node voltage angle in 3 x nnode matrix
% Vopt - node voltage phasor in 3 x nnode matrix
% Sopt - line complex power in 3 x nline matrix
% wopt - node DER complex dispatch in 3 x nnode matrix
% demopt - node complex loads in 3 x nnode matrix
% sopt - node total complex power in 3 x nnode matrix

base = network.base;
nodes = network.nodes;
lines = network.lines;
configs = network.configs;
loads = network.loads;
caps = network.caps;
cons = network.cons;

nnode = nodes.nnode;
nline = lines.nline;

% Separate portions of X corresponding to each phase
Xa = X(1:nvar);
Xb = X(nvar+1:2*nvar);
Xc = X(2*nvar+1:3*nvar);

% Place in 3 x nvar matrix
XX = [Xa Xb Xc]';

% Portion of XX corresponding to the squared voltage magnitude
Eopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

% Portion of XX corresponding to the voltage angle
Dopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

% Compute optimal voltage phasor from Eopt, Dopt
Vopt = sqrt(Eopt).*exp(1j*pi/180*Dopt).*nodes.PH;

% Compute optimal line current from voltage phasors
for k1 = 1:lines.nline
    Iopt(:,k1) = lines.FYpu(:,:,k1)*(Vopt(:,lines.TXnum(k1)) - Vopt(:,lines.RXnum(k1)));
end
Iopt = Iopt.*lines.PH;

% Portion of XX corresponding to line real power
Popt = XX(:,1:nline);
XX(:,1:nline) = [];

% Portion of XX corresponding to line reactive power
Qopt = XX(:,1:nline);
XX(:,1:nline) = [];

% Line complex power
Sopt = Popt + 1j*Qopt;

% Portion of XX corresponding to node DER real disptach
uopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

% Portion of XX corresponding to node DER reactive disptach
vopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

% Node complex DER dispatch
wopt = uopt + 1j*vopt;

% Node loads (demand)
demopt = loads.spu.*(loads.aPQ + loads.aZ.*Eopt).*nodes.PH;

% Node total power
for k1 = 1:nnode
    sopt(:,k1) = sum(Sopt(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
        sum(Sopt(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2);
end
% sopt = dempot + wopt - 1j*caps.cappu;

end