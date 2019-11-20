function [nvar, Aineq, bineq, Aeq, beq] = create_linear_constraints_taylor_expansion(network,slacknode,Vslack,sim)
%%

% Michael Sankur - msankur@lbl.gov
% 2018.06.01

% Thius function formualtes the equations of the linearized unbalanced
% power flow model into matrix/vector format where:
% Aineq X <= bineq
% and
% Aeq X = beq
% so these can be used as inequality and equality constraints for
% optimization programs

% INPUT(S)
% network - struct containing all pertinent network information,
% including all other structs
% base - struct containing base values
% nodes - struct containing node parameters
% lines - struct containing line parameters
% loads - struct containing load parameters
% caps - struct containing capacitor parameters
% cons - struct containing controller parameters
% slacknode - index of slack node
% Vslack - voltage reference for slack node
% sim - struct containing parameters pertinent to the simulation, such as
% objective function coefficients, higher order terms obtained by solving
% power flow

% OUTPUT(S)
% nvar - number of variables per phase
% Aineq - matrix for all inequality constraints Aineq X <= bineq
% bineq - RHS vector for all inequality constrains Aineq X <= bineq
% Aeq - matrix for all equality constraints Aeq X = beq
% beq - RHS vector for all equality constrains Aeq X = beq

% The variable names are as follows:
% E - squared voltage magnitude where E_n^phi = | V_n^phi |^2 [pu]^2
% T - voltage angle [deg]
% P - line real power [pu]
% Q - line reactive power [pu]
% u - controllable node real power [pu]
% v - controllable node reactive power [pu]

% The variable X in the inequality and equality constraints is
% X = [X_a X_b X_c]
% with
% X_phi = [E^phi T^phi P^phi Q^phi u^phi v^phi]
% with
% E^phi = [E_1^phi E_2^phi ... E_nnode^phi]
% T^phi = [T_1^phi T_2^phi ... T_nnode^phi]
% P^phi = [P_1^phi P_2^phi ... P_nline^phi]
% Q^phi = [Q_1^phi Q_2^phi ... Q_nline^phi]
% u^phi = [u_1^phi u_2^phi ... u_nnode^phi]
% v^phi = [v_1^phi v_2^phi ... v_nnode^phi]


%%

base = network.base;
nodes = network.nodes;
lines = network.lines;
configs = network.configs;
loads = network.loads;
caps = network.caps;
cons = network.cons;

% node parameters
nnode = nodes.nnode;
NPH = nodes.PH;
inlines = nodes.inlines;
innodes = nodes.innodes;
outlines = nodes.outlines;
outnodes = nodes.outnodes;

% line paramters
LPH = lines.PH;
nline = lines.nline;
TXnum = lines.TXnum;
RXnum = lines.RXnum;
FZpu = lines.FZpu;

% load parameters
spu = loads.spu;
aS = loads.aPQ;
aI = loads.aI;
aZ = loads.aZ;

% capacitor parameters
cappu = caps.cappu;

% controller parameters
wmaxpu = cons.wmaxpu;

% sim parameters
% rho = sim.rho;

% iteration parameters
if isempty(sim.VNR)
    VNR = [1*ones(1,nnode);
        1*exp(1j*-120*pi/180)*ones(1,nnode);
        1*exp(1j*120*pi/180)*ones(1,nnode)].*NPH;
else
    VNR = sim.VNR;
end
% L_mn
if isempty(sim.Lmn)
    Lmn = zeros(3,nline);
else
    Lmn = sim.Lmn;
end
% H_mn
if isempty(sim.Hmn)
    Hmn = zeros(3,nline);
else
    Hmn = sim.Hmn;
end

% number of variables per phase
nvar = nnode + nnode + nline + nline + nnode + nnode;

%%

Abar = real(sim.Vm(1,:));
Bbar = imag(sim.Vm(1,:));

Cbar = real(sim.Imn(1,:));
Dbar = imag(sim.Imn(1,:));

ubar = real(sim.w(1,:));
vbar = imag(sim.w(1,:));

%% Set up optimization matrices

Aineq = []; % Inequality constraint matrix
bineq = []; % Inequality constraint vector
Aeq = []; % Equality constraint matrix
beq = []; % Equality constraint vector

%% Network slack node voltage magnitude constraints

tempA = zeros(1,nnode); tempA(1,slacknode) = 2*Abar(slacknode);
tempB = zeros(1,nnode); tempB(1,slacknode) = 2*Bbar(slacknode);
tempCmn = zeros(1,nline);
tempDmn = zeros(1,nline);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

Aeq = [Aeq;
    tempA tempB tempCmn tempDmn tempu tempv];

beq = [beq;
    Vslack(1)^2 - Abar(slacknode)^2 - Bbar(slacknode)^2];

%% Network voltage magnitude constraints

for kn = 1:nnode
    
    tempA = zeros(1,nnode); tempA(1,kn) = -2*Abar(kn);
    tempB = zeros(1,nnode); tempB(1,kn) = -2*Bbar(kn);
    tempCmn = zeros(1,nline);
    tempDmn = zeros(1,nline);
    tempu = zeros(1,nnode);
    tempv = zeros(1,nnode);
    
    Aineq = [Aineq;
        tempA tempB tempCmn tempDmn tempu tempv];
    
    bineq = [bineq;
        Abar(kn)^2 + Bbar(kn)^2 - 0.95^2];
    
    
    tempA = zeros(1,nnode); tempA(1,kn) = 2*Abar(kn);
    tempB = zeros(1,nnode); tempB(1,kn) = 2*Bbar(kn);
    tempCmn = zeros(1,nline);
    tempDmn = zeros(1,nline);
    tempu = zeros(1,nnode);
    tempv = zeros(1,nnode);
    
    Aineq = [Aineq;
        tempA tempB tempCmn tempDmn tempu tempv];
    
    bineq = [bineq;
         1.05^2 - Abar(kn)^2 - Bbar(kn)^2];    
    
end

%% Power Flow - No Losses Explicity Modeled

for kn = 1:nnode

    % real power
    tempA = zeros(1,nnode);
    tempA(1,kn) = -sum(Cbar(inlines(inlines(:,kn) ~= 0,kn))) ...
        + sum(Cbar(outlines(outlines(:,kn) ~= 0,kn))) ...
        + real(spu(1,kn))*(aI(1,kn)*(Abar(kn)^2 + Bbar(kn)^2)^(-1/2)*Abar(kn)) ...
        + real(spu(1,kn))*(aZ(1,kn)*2*Abar(kn));
    tempB = zeros(1,nnode);
    tempB(1,kn) = -sum(Dbar(inlines(inlines(:,kn) ~= 0,kn))) ...
        + sum(Dbar(outlines(outlines(:,kn) ~= 0,kn))) ...
        + real(spu(1,kn))*(aI(1,kn)*(Abar(kn)^2 + Bbar(kn)^2)^(-1/2)*Bbar(kn)) ...
        + real(spu(1,kn))*(aZ(1,kn)*2*Bbar(kn));
    tempCmn = zeros(1,nline);
    tempCmn(1,inlines(inlines(:,kn)~=0,kn)) = -Abar(kn);
    tempCmn(1,outlines(outlines(:,kn)~=0,kn)) = Abar(kn);
    tempDmn = zeros(1,nline);
    tempDmn(1,inlines(inlines(:,kn)~=0,kn)) = -Bbar(kn);
    tempDmn(1,outlines(outlines(:,kn)~=0,kn)) = Bbar(kn);
    tempu = zeros(1,nnode);
    tempu(kn) = 1;
    tempv = zeros(1,nnode);
    
    Aeq = [Aeq;
        tempA tempB tempCmn tempDmn tempu tempv];
    
    beq = [beq;
        sum(Abar(kn)*Cbar(inlines(inlines(:,kn) ~= 0,kn)) + Bbar(kn)*Dbar(inlines(inlines(:,kn) ~= 0,kn))) ...
        - sum(Abar(kn)*Cbar(outlines(outlines(:,kn) ~= 0,kn)) + Bbar(kn)*Dbar(outlines(outlines(:,kn) ~= 0,kn))) ...
        - real(spu(1,kn))*(aS(1,kn) + aI(1,kn)*(Abar(kn)^2 + Bbar(kn)^2)^(1/2) + aZ(1,kn)*(Abar(kn)^2 + Bbar(kn)^2)) - ubar(kn)];
    
    % reactive power
    tempA = zeros(1,nnode);
    tempA(1,kn) = sum(Dbar(inlines(inlines(:,kn) ~= 0,kn))) ...
        - sum(Dbar(outlines(outlines(:,kn) ~= 0,kn))) ...
        + imag(spu(1,kn))*(aI(1,kn)*(Abar(kn)^2 + Bbar(kn)^2)^(-1/2)*Abar(kn)) ...
        + imag(spu(1,kn))*(aZ(1,kn)*2*Abar(kn));
    tempB = zeros(1,nnode);
    tempB(1,kn) = -sum(Cbar(inlines(inlines(:,kn) ~= 0,kn))) ...
        + sum(Cbar(outlines(outlines(:,kn) ~= 0,kn))) ...
        + imag(spu(1,kn))*(aI(1,kn)*(Abar(kn)^2 + Bbar(kn)^2)^(-1/2)*Bbar(kn)) ...
        + imag(spu(1,kn))*(aZ(1,kn)*2*Bbar(kn));
    tempCmn = zeros(1,nline);
    tempCmn(1,inlines(inlines(:,kn)~=0,kn)) = -Bbar(kn);
    tempCmn(1,outlines(outlines(:,kn)~=0,kn)) = Bbar(kn);
    tempDmn = zeros(1,nline);
    tempDmn(1,inlines(inlines(:,kn)~=0,kn)) = Abar(kn);
    tempDmn(1,outlines(outlines(:,kn)~=0,kn)) = -Abar(kn);
    tempu = zeros(1,nnode);
    tempv = zeros(1,nnode);
    tempv(kn) = 1;
    
    Aeq = [Aeq;
        tempA tempB tempCmn tempDmn tempu tempv];
    
    beq = [beq;
        sum(-Abar(kn)*Dbar(inlines(inlines(:,kn) ~= 0,kn)) + Bbar(kn)*Cbar(inlines(inlines(:,kn) ~= 0,kn))) ...
        - sum(-Abar(kn)*Dbar(outlines(outlines(:,kn) ~= 0,kn)) + Bbar(kn)*Cbar(outlines(outlines(:,kn) ~= 0,kn))) ...
        - imag(spu(1,kn))*(aS(1,kn) + aI(1,kn)*(Abar(kn)^2 + Bbar(kn)^2)^(1/2) + aZ(1,kn)*(Abar(kn)^2 + Bbar(kn)^2)) - vbar(kn) + cappu(1,kn)];

end


%%

for kn = 1:nnode
    
    if wmaxpu(1,kn) == 0
        tempA = zeros(1,nnode);
        tempB = zeros(1,nnode);
        tempCmn = zeros(1,nline);
        tempDmn = zeros(1,nline);
        tempu = zeros(1,nnode); tempu(kn) = 1;
        tempv = zeros(1,nnode);

        Aeq = [Aeq;
            tempA tempB tempCmn tempDmn tempu tempv];

        beq = [beq;
            0];

        tempA = zeros(1,nnode);
        tempB = zeros(1,nnode);
        tempCmn = zeros(1,nline);
        tempDmn = zeros(1,nline);
        tempu = zeros(1,nnode);
        tempv = zeros(1,nnode); tempv(kn) = 1;

        Aeq = [Aeq;
            tempA tempB tempCmn tempDmn tempu tempv];

        beq = [beq;
            0];
    end

    if wmaxpu(1,kn) > 0
        
        tempA = zeros(1,nnode);
        tempB = zeros(1,nnode);
        tempCmn = zeros(1,nline);
        tempDmn = zeros(1,nline);
        tempu = zeros(1,nnode); tempu(kn) = 2*ubar(kn);
        tempv = zeros(1,nnode); tempv(kn) = 2*vbar(kn);

        Aineq = [Aineq;
            tempA tempB tempCmn tempDmn tempu tempv];

        bineq = [bineq;
            wmaxpu(1,kn)^2 - ubar(kn)^2 - vbar(kn)^2];

    end
    
    
end

end