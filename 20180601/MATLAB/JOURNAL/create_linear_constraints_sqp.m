function [nvar, Aineq, bineq, Aeq, beq] = create_linear_constraints_sqp(network,slacknode,Vslack,sim)
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
% X = [X_a; X_b; X_c]
% with
% X_phi = [E^phi; T^phi; P^phi; Q^phi; u^phi; v^phi]
% with
% E^phi = [E_1^phi E_2^phi ... E_nnode^phi]^T
% T^phi = [T_1^phi T_2^phi ... T_nnode^phi]^T
% P^phi = [P_1^phi P_2^phi ... P_nline^phi]^T
% Q^phi = [Q_1^phi Q_2^phi ... Q_nline^phi]^T
% u^phi = [u_1^phi u_2^phi ... u_nnode^phi]^T
% v^phi = [v_1^phi v_2^phi ... v_nnode^phi]^T


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
aPQ = loads.aPQ;
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

refnode = 'TX';

%% Set up optimization matrices

Aineq = []; % Inequality constraint matrix
bineq = []; % Inequality constraint vector
Aeq = []; % Equality constraint matrix
beq = []; % Equality constraint vector

%% Network slack node voltage magnitude constraints

% Set slacknode square voltage magnitude

tempE = zeros(1,nnode); tempE(1,slacknode) = 1;
tempT = zeros(1,nnode);
tempP = zeros(1,nline);
tempQ = zeros(1,nline);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

Aeq = [Aeq;
    tempE tempT tempP tempQ tempu tempv, zeros(1,nvar), zeros(1,nvar);
    zeros(1,nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,nvar);
    zeros(1,nvar), zeros(1,nvar), tempE tempT tempP tempQ tempu tempv];

beq = [beq;
    abs(Vslack)];

clear tempE tempT tempP tempQ tempu tempv

%% Network slack node voltage angle constraints

% Set slacknode square voltage angle

tempE = zeros(1,nnode);
tempT = zeros(1,nnode); tempT(1,slacknode) = 1;
tempP = zeros(1,nline);
tempQ = zeros(1,nline);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

Aeq = [Aeq;
    tempE tempT tempP tempQ tempu tempv, zeros(1,nvar), zeros(1,nvar);
    zeros(1,nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,nvar);
    zeros(1,nvar), zeros(1,nvar), tempE tempT tempP tempQ tempu tempv];

beq = [beq;
    180/pi*angle(Vslack)];

clear tempE tempT tempP tempQ tempu tempv

%% Voltage magnitude equality and inequality constraints

% If no phase present at node - Set voltage magnitude to zero
% If phase present at node - Add inequality constraint on squared voltage
% magnitude - 0.95^2 <= E_n^phi <= 1.05^2

for ph = 1:3
    for k1 = 2:nnode
% If no phase present at node - Set voltage magnitude to zero        
        if NPH(ph,k1) == 0
            tempE = zeros(1,nnode); tempE(1,k1) = 1;
            tempT = zeros(1,nnode);
            tempP = zeros(1,nline);
            tempQ = zeros(1,nline);
            tempu = zeros(1,nnode);
            tempv = zeros(1,nnode);
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0];
% If phase present at node - Add inequality constraint on squared voltage
% magnitude - 0.95^2 <= E_n^phi <= 1.05^2
        elseif NPH(ph,k1) == 1
            tempE = zeros(2,nnode); tempE(:,k1) = [1; -1];
            tempT = zeros(2,nnode);
            tempP = zeros(2,nline);
            tempQ = zeros(2,nline);
            tempu = zeros(2,nnode);
            tempv = zeros(2,nnode);
            Aineq = [Aineq;
                zeros(2,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(2,(3-ph)*nvar)];
            bineq = [bineq;
                1.05^2;
                -(0.95^2)];
        end
    end
end

clear tempE tempT tempP tempQ tempu tempv

%% Approximate |V| equality constraints

% |V_n^phi| = sqrt(E_n^phi)
% Approximate |V_n^phi| using first order Taylor Expansion
% |V_n^phi| ~ 0.5*(1 + E_n^phi)
% for ph = 1:3
%     for k1 = 1:n
%         if PH(ph,k1) == 0
%             tempVmag = zeros(1,n); tempVmag(k1) = 1;
%             Aeq = [Aeq;
%                 zeros(1,(ph-1)*nvar), zeros(1,n) tempVmag zeros(1,5*n), zeros(1,(3-ph)*nvar)];
%             beq = [beq;
%                 0];
%         elseif PH(ph,k1) == 1
%             tempE = zeros(1,n); tempE(k1) = -0.5;
%             tempVmag = zeros(1,n); tempVmag(k1) = 1;%             
%             Aeq = [Aeq;
%                 zeros(1,(ph-1)*nvar), tempE tempVmag zeros(1,5*n), zeros(1,(3-ph)*nvar)];
%             beq = [beq;
%                 0.5];
%         end
%     end
% end

%% Power flow equality constraints

% FIXXXXX COMMENTS

% If no phase present at node - Set power entering nodes to zero -
% P_m^phi = 0, Q_m^phi = 0
% If phase present at node - 
% Real
% sum_l P_lm^phi - u_m^phi - sum_n ( P_mn^phi ) = p_n^phi + sum_n Re{L_mn^phi}
% p_n^phi = p_n^phi*(A_PQ,n^phi + A_Z,n^phi E_n^phi) + u_n^phi
% Imag
% sum_l Q_lm^phi - v_m^phi - sum_n ( Q_mn^phi ) = q_n^phi - c_n^phi + sum_n Im{L_mn^phi}
% q_n^phi = q_n^phi*(A_PQ,n^phi + A_Z,n^phi E_n^phi) - cap_n^phi + v_n^phi

tempE = zeros(1,nnode);
tempT = zeros(1,nnode);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

% If no phase phi present on line (m,n) - Set power on phase to zero:
% P_mn^phi = 0, Q_mn^phi = 0
for ph = 1:3
    for k1 = 1:nline
        if LPH(ph,k1) == 0
            tempPQ = zeros(1,nline); tempPQ(1,k1) = 1;
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), tempE tempT tempPQ zeros(1,nline) tempu tempv, zeros(1,(3-ph)*nvar);
                zeros(1,(ph-1)*nvar), tempE tempT zeros(1,nline) tempPQ tempu tempv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0;
                0];
        end
    end
end

for ph = 1:3
    for k1 = 2:nnode
% If no phase phi present at node m - Set power entering and leaving nodes
% on phase to zero:
% P_lm^phi = 0, Q_lm^phi = 0, P_mn^phi = 0, Q_mn^phi = 0
        if NPH(ph,k1) == 0
            for k2 = 1:size(inlines,1)
                if inlines(k2,k1) ~= 0
                    tempE = zeros(1,nnode);
                    tempT = zeros(1,nnode);
                    tempPQ = zeros(1,nline); tempPQ(1,inlines(k2,k1)) = 1;
                    tempu = zeros(1,nnode);
                    tempv = zeros(1,nnode);
                    Aeq = [Aeq;
                        zeros(1,(ph-1)*nvar), tempE tempT tempPQ zeros(1,nline) tempu tempv, zeros(1,(3-ph)*nvar);
                        zeros(1,(ph-1)*nvar), tempE tempT zeros(1,nline) tempPQ tempu tempv, zeros(1,(3-ph)*nvar)];
                    beq = [beq;
                        0;
                        0];
                end
            end
            for k2 = 1:size(outlines,1)
                if outlines(k2,k1) ~= 0
                    tempE = zeros(1,nnode);
                    tempT = zeros(1,nnode);
                    tempPQ = zeros(1,nline); tempPQ(1,outlines(k2,k1)) = 1;
                    tempu = zeros(1,nnode);
                    tempv = zeros(1,nnode);
                    Aeq = [Aeq;
                            zeros(1,(ph-1)*nvar), tempE tempT tempPQ zeros(1,nline) tempu tempv, zeros(1,(3-ph)*nvar);
                            zeros(1,(ph-1)*nvar), tempE tempT zeros(1,nline) tempPQ tempu tempv, zeros(1,(3-ph)*nvar)];
                    beq = [beq;
                        0;
                        0];
                end
            end
% If phase present at node - Conservation of power with ZIP demand model
% Real
% sum_(l:(l,m) in Edges) P_lm^phi - sum_n ( P_mn^phi ) - p_n^phi = 0
% p_n^phi = Re{d_n^phi}*(A_PQ,n^phi + A_I,n^phi*(1 + E_n^phi)/2 + A_Z,n^phi*E_n^phi) + u_n^phi
% Imag
% sum_l Q_lm^phi - v_m^phi - sum_n ( Q_mn^phi ) - q_n^phi = - c_n^phi
% q_n^phi = Im{d_n^phi}*(A_PQ,n^phi + A_I,n^phi*(1 + E_n^phi)/2 + A_Z,n^phi*E_n^phi) + v_n^phi
        elseif NPH(ph,k1) == 1
            tempE = zeros(1,nnode); tempE(1,k1) = -spu(ph,k1)*(aI(ph,k1)/2 + aZ(ph,k1));
            tempT = zeros(1,nnode);
            tempPQ = zeros(1,nline);
            tempPQ(inlines(inlines(:,k1) ~= 0,k1)) = 1;
            tempPQ(outlines(outlines(:,k1) ~= 0,k1)) = -1;
            tempuv = zeros(1,nnode); tempuv(1,k1) = -1;
            Aeq  = [Aeq;
                zeros(1,(ph-1)*nvar), real(tempE) tempT tempPQ zeros(1,nline) tempuv zeros(1,nnode), zeros(1,(3-ph)*nvar);
                zeros(1,(ph-1)*nvar), imag(tempE) tempT zeros(1,1*nline) tempPQ zeros(1,nnode) tempuv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                real(spu(ph,k1))*(aPQ(ph,k1) + aI(ph,k1)/2);
                imag(spu(ph,k1))*(aPQ(ph,k1) + aI(ph,k1)/2) - cappu(ph,k1)];
            if strcmp(refnode,'TX') == 1
                beq(end-1:end) = beq(end-1:end) + ...
                    [sum(real(Lmn(ph,inlines(inlines(:,k1)~=0,k1))),2);
                    sum(imag(Lmn(ph,inlines(inlines(:,k1)~=0,k1))),2)];
            end
            if strcmp(refnode,'RX') == 1
                beq(end-1:end) = beq(end-1:end) + ...
                    [sum(real(Lmn(ph,outlines(outlines(:,k1)~=0,k1))),2);
                    sum(imag(Lmn(ph,outlines(outlines(:,k1)~=0,k1))),2)];
            end
%             outlines(outlines(:,k1) ~= 0,k1);
%             beq(end-1) = beq(end-1) + real(sum(Lmn(ph,outlines(outlines(ph,k1) ~= 0,k1)),2));
%             beq(end) = beq(end) + imag(sum(Lmn(ph,outlines(outlines(ph,k1) ~= 0,k1)),2));
        end        
    end        
end

clear tempE tempT tempP tempQ tempPQ tempu tempv tempuv

% sum(sim.Lmn(:,outlines(outlines(:,12)~=0,12)),2)

%% Calculate gamma terms for each node

gammaFZpu = zeros(3,3,nline);
for k1 = 1:nline
%     txnode = TXnum(k1);
%     rxnode = RXnum(k1);
%     node = rxnode;
    if strcmp(refnode,'TX') == 1
        node = TXnum(k1);
    end
    if strcmp(refnode,'RX') == 1
        node = RXnum(k1);
    end
    gamman = VNR(:,node)*(1./VNR(:,node)).';
    if NPH(1,node) == 0
        gamman(1,:) = 0; gamman(:,1) = 0;
    end
    if NPH(2,node) == 0
        gamman(2,:) = 0; gamman(:,2) = 0;
    end
    if NPH(3,node) == 0
        gamman(3,:) = 0; gamman(:,3) = 0;
    end
    gamman;
    gammaFZpu(:,:,k1) = gamman.*conj(FZpu(:,:,k1));
end
Mmn = real(gammaFZpu);
Nmn = imag(gammaFZpu);

%% Voltage magnitude approximation equality constraints

% If no phase present at node - Set voltage at phase to zero - E_n^phi = 0
% If phase present at node - Propagate voltage up feeder (toward head)
% -E_m + E_n + 2 M_mn P_mn + 2 N_mn Q_mn = -H_mn

tempT = zeros(1,nnode);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

% Phase A
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if NPH(1,txnode) == 1 && NPH(1,rxnode) == 1 && LPH(1,k1) == 1
        tempE = zeros(1,nnode); tempE(1,[txnode rxnode]) = [-1 1];
        tempPa = zeros(1,nline); tempPa(1,k1) = 2*Mmn(1,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = 2*Mmn(1,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = 2*Mmn(1,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = -2*Nmn(1,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = -2*Nmn(1,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = -2*Nmn(1,3,k1);

        Aeq = [Aeq;
            tempE tempT tempPa tempQa tempu tempv, ...
            zeros(1,nnode) tempT tempPb tempQb tempu tempv, ...
            zeros(1,nnode) tempT tempPc tempQc tempu tempv];
        beq = [beq;
            0];
        if strcmp(refnode,'TX') == 1
            beq(end) = beq(end) + Hmn(1,k1);
        end
        if strcmp(refnode,'RX') == 1
            beq(end) = beq(end) - Hmn(1,k1);
        end
    end
end

% Phase B
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);
    
    if NPH(2,txnode) == 1 && NPH(2,rxnode) == 1 && LPH(2,k1) == 1
        tempE = zeros(1,nnode); tempE(1,[txnode rxnode]) = [-1 1];
        tempPa = zeros(1,nline); tempPa(1,k1) = 2*Mmn(2,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = 2*Mmn(2,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = 2*Mmn(2,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = -2*Nmn(2,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = -2*Nmn(2,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = -2*Nmn(2,3,k1);

        Aeq = [Aeq;
            zeros(1,nnode) tempT tempPa tempQa tempu tempv, ...
            tempE tempT tempPb tempQb tempu tempv, ...
            zeros(1,nnode) tempT tempPc tempQc tempu tempv];
        beq = [beq;
            0];
        if strcmp(refnode,'TX') == 1
            beq(end) = beq(end) + Hmn(2,k1);
        end
        if strcmp(refnode,'RX') == 1
            beq(end) = beq(end) - Hmn(2,k1);
        end
    end
end

% Phase C
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if NPH(3,txnode) == 1 && NPH(3,rxnode) == 1 && LPH(3,k1) == 1
        tempE = zeros(1,nnode); tempE(1,[txnode rxnode]) = [-1 1];
        tempPa = zeros(1,nline); tempPa(1,k1) = 2*Mmn(3,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = 2*Mmn(3,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = 2*Mmn(3,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = -2*Nmn(3,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = -2*Nmn(3,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = -2*Nmn(3,3,k1);

        Aeq = [Aeq;
            zeros(1,nnode) tempT tempPa tempQa tempu tempv, ...
            zeros(1,nnode) tempT tempPb tempQb tempu tempv, ...
            tempE tempT tempPc tempQc tempu tempv];
        beq = [beq;
            0];
        if strcmp(refnode,'TX') == 1
            beq(end) = beq(end) + Hmn(3,k1);
        end
        if strcmp(refnode,'RX') == 1
            beq(end) = beq(end) - Hmn(3,k1);
        end
    end
end

clear tempE tempT tempu tempv
clear tempPa tempQa tempPb tempQb tempPc tempQc

%% Voltage angle approximation equality constraints

% If no phase present at node - No constraint on phase angle
% If phase present at node - Propagate voltage up feeder (toward head)
% |V_m||V_n|(-Theta_m + Theta_n) - N_mn P_mn - M_mn Q_mn = 0
% Above is not exact equation, read papers to find it

tempE = zeros(1,nnode);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

% Phase A
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if NPH(1,txnode) == 1 && NPH(1,rxnode) == 1 && LPH(1,k1) == 1
        tempT = zeros(1,nnode); tempT(1,[txnode rxnode]) = pi/180*[-1 1];
%         magprod = abs(VNR(1,txnode))*abs(VNR(1,rxnode))
        tempT = tempT*abs(VNR(1,txnode))*abs(VNR(1,rxnode));
%         angTX = 180/pi*angle(VNR(1,txnode))
%         angRX = 180/pi*angle(VNR(1,rxnode))
%         cosangdif = cos(angle(VNR(1,txnode)) - angle(VNR(1,rxnode)))
%         tempT = tempT*cos(angle(VNR(1,txnode)) - angle(VNR(1,rxnode)));
        tempPa = zeros(1,nline); tempPa(1,k1) = -Nmn(1,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = -Nmn(1,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = -Nmn(1,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = -Mmn(1,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = -Mmn(1,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = -Mmn(1,3,k1);

        Aeq = [Aeq;
            tempE tempT tempPa tempQa tempu tempv, ...
            tempE zeros(1,nnode) tempPb tempQb tempu tempv, ...
            tempE zeros(1,nnode) tempPc tempQc tempu tempv];
        beq = [beq;
            0];
%         rhs = abs(VNR(1,txnode))*abs(VNR(1,rxnode))*sin(angle(VNR(1,txnode)) - angle(VNR(1,rxnode)))
%         beq(end) = beq(end) + abs(VNR(1,txnode))*abs(VNR(1,rxnode))*sin(angle(VNR(1,txnode)) - angle(VNR(1,rxnode)));
    end
end

% Phase B
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if NPH(2,txnode) == 1 && NPH(2,rxnode) == 1 && LPH(2,k1) == 1
        tempT = zeros(1,nnode); tempT(1,[txnode rxnode]) = pi/180*[-1 1];
        tempT = tempT*abs(VNR(2,txnode))*abs(VNR(2,rxnode));
%         cos(angle(VNR(2,txnode)) - angle(VNR(2,rxnode)))
%         tempT = tempT*cos(angle(VNR(2,txnode)) - angle(VNR(2,rxnode)));
        tempPa = zeros(1,nline); tempPa(1,k1) = -Nmn(2,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = -Nmn(2,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = -Nmn(2,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = -Mmn(2,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = -Mmn(2,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = -Mmn(2,3,k1);

        Aeq = [Aeq;
            tempE zeros(1,nnode) tempPa tempQa tempu tempv, ...
            tempE tempT tempPb tempQb tempu tempv, ...
            tempE zeros(1,nnode) tempPc tempQc tempu tempv];
        beq = [beq;
            0];
%         abs(VNR(2,txnode))*abs(VNR(2,rxnode))*sin(angle(VNR(2,txnode)) - angle(VNR(2,rxnode)))
%         beq(end) = beq(end) + abs(VNR(2,txnode))*abs(VNR(2,rxnode))*sin(angle(VNR(2,txnode)) - angle(VNR(2,rxnode)));
    end
end

% Phase C
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if NPH(3,txnode) == 1 && NPH(3,rxnode) == 1 && LPH(3,k1) == 1
        tempT = zeros(1,nnode); tempT(1,[txnode rxnode]) = pi/180*[-1 1];
        tempT = tempT*abs(VNR(3,txnode))*abs(VNR(3,rxnode));
%         cos(angle(VNR(3,txnode)) - angle(VNR(3,rxnode)))
%         tempT = tempT*cos(angle(VNR(3,txnode)) - angle(VNR(3,rxnode)));
        tempPa = zeros(1,nline); tempPa(1,k1) = -Nmn(3,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = -Nmn(3,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = -Nmn(3,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = -Mmn(3,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = -Mmn(3,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = -Mmn(3,3,k1);

        Aeq = [Aeq;
            tempE zeros(1,nnode) tempPa tempQa tempu tempv, ...
            tempE zeros(1,nnode) tempPb tempQb tempu tempv, ...
            tempE tempT tempPc tempQc tempu tempv];
        beq = [beq;
            0];
%         abs(VNR(3,txnode))*abs(VNR(3,rxnode))*sin(angle(VNR(3,txnode)) - angle(VNR(3,rxnode)))
%         beq(end) = beq(end) + abs(VNR(3,txnode))*abs(VNR(3,rxnode))*sin(angle(VNR(3,txnode)) - angle(VNR(3,rxnode)));
    end
end

clear tempT tempPa tempQa tempPb tempQb tempPc tempQc

%% Zero control outputs for nonexistent phases

% For all phases at all nodes that are non control, set inverter output to
% zero - u_n^phi = 0, v_n^phi = 0

tempE = zeros(1,nnode);
tempT = zeros(1,nnode);
tempP = zeros(1,nline);
tempQ = zeros(1,nline);

for ph = 1:3
    for k1 = 1:nnode
        if NPH(ph,k1) == 0
            tempu = zeros(1,nnode); tempu(1,k1) = 1;
            tempv = zeros(1,nnode);
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0];
            tempu = zeros(1,nnode); 
            tempv = zeros(1,nnode); tempv(1,k1) = 1;
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0];
        end
    end
end

clear tempE tempT tempP tempQ tempu tempu tempv

%% Zero control real power dispatch

% For cons with zero real power dispatch, set output to zero
% u_n^phi = 0

tempE = zeros(1,nnode);
tempT = zeros(1,nnode);
tempP = zeros(1,nline);
tempQ = zeros(1,nline);

if isfield(cons,'realpower') == 1
    for ph = 1:3
        for k1 = 1:nnode
            if cons.realpower(ph,k1) == 0
                tempu = zeros(1,nnode); tempu(1,k1) = 1;
                tempv = zeros(1,nnode);
                Aeq = [Aeq;
                    zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
                beq = [beq;
                    0];
            elseif cons.realpower(ph,k1) == 1

            end
        end
    end
end

clear tempE tempT tempP tempQ tempu tempu tempv

%% Zero control reactive power dispatch

% For cons with zero reactive power dispatch, set output to zero
% u_n^phi = 0

tempE = zeros(1,nnode);
tempT = zeros(1,nnode);
tempP = zeros(1,nline);
tempQ = zeros(1,nline);

if isfield(cons,'reactivepower') == 1
    for ph = 1:3
        for k1 = 1:nnode
            if cons.reactivepower(ph,k1) == 0
                tempu = zeros(1,nnode);
                tempv = zeros(1,nnode); tempv(1,k1) = 1;
                Aeq = [Aeq;
                    zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
                beq = [beq;
                    0];            
            elseif cons.reactivepower(ph,k1) == 1

            end
        end
    end
end

clear tempE tempT tempP tempQ tempu tempu tempv

%% Control negative real power dispatch (only real power generation)

% For cons with negative real power dispatch, set output to
% negative. This represents a controller only outputtting powre
% (generation)
% u_n^phi <= 0 

tempE = zeros(1,nnode);
tempT = zeros(1,nnode);
tempP = zeros(1,nline);
tempQ = zeros(1,nline);

if isfield(cons,'realgen') == 1
    for ph = 1:3
        for k1 = 1:nnode
            if cons.realgen(ph,k1) == 1
                tempu = zeros(1,nnode); tempu(1,k1) = 1;
                tempv = zeros(1,nnode);
                Aineq = [Aineq;
                    zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
                bineq = [bineq;
                    0];
            end
        end
    end
end

clear tempE tempT tempP tempQ tempu tempu tempv

%% Control mimimum and/or maximum real power dispatch

% For cons with bounds on the minimum and/or maximum real power
% dispatch. If applicable:
% u_n^phi <= u_n^phi(max)
% u_n^phi(max) <= u_n^phi

tempE = zeros(1,nnode);
tempT = zeros(1,nnode);
tempP = zeros(1,nline);
tempQ = zeros(1,nline);

if isfield(cons,'asdf') == 1
    for ph = 1:3
        for k1 = 1:nnode
            if cons.negativerealpower(ph,k1) == 1
                tempu = zeros(1,nnode); tempu(1,k1) = -1;
                tempv = zeros(1,nnode);
                Aineq = [Aineq;
                    zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
                bineq = [bineq;
                    uminpu(ph,k1)];
            end
            if cons.negativerealpower(ph,k1) == 1
                tempu = zeros(1,nnode); tempu(1,k1) = 1;
                tempv = zeros(1,nnode);
                Aineq = [Aineq;
                    zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
                bineq = [bineq;
                    umaxpu(ph,k1)];
            end
        end
    end
end

clear tempE tempT tempP tempQ tempu tempu tempv

%% Control mimimum and/or maximum reactive power dispatch

% For cons with bounds on the minimum and/or maximum reactive power
% dispatch. If applicable:
% v_n^phi <= v_n^phi(max)
% v_n^phi(max) <= v_n^phi

tempE = zeros(1,nnode);
tempT = zeros(1,nnode);
tempP = zeros(1,nline);
tempQ = zeros(1,nline);

if isfield(cons,'asdf') == 1
    for ph = 1:3
        for k1 = 1:nnode
            if cons.negativerealpower(ph,k1) == 1
                tempu = zeros(1,nnode); 
                tempv = zeros(1,nnode); tempv(1,k1) = -1;
                Aineq = [Aineq;
                    zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
                bineq = [bineq;
                    vminpu(ph,k1)];
            end
            if cons.negativerealpower(ph,k1) == 1
                tempu = zeros(1,nnode); 
                tempv = zeros(1,nnode); tempv(1,k1) = 1;
                Aineq = [Aineq;
                    zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
                bineq = [bineq;
                    vmaxpu(ph,k1)];
            end
        end
    end
end

clear tempE tempT tempP tempQ tempu tempu tempv

%% Control apparent power constraint - circle

% Impose constraint on apparent power using linear approximation of L2
% norm.
% cos(delta) u_n^phi + sin(delta) v_n^phi <= wmaxpu_n^phi
% delta = 0:dx:2*pi

% ncirc = 12;
% delta = linspace(0,2*pi,ncirc+1);
% 
% if isfield(cons,'wmaxpu') == 1
%     for ph = 1:3
%         for k1 = 1:nnode
%             if wmaxpu(ph,k1) == 0
%                 tempE = zeros(1,nnode);
%                 tempT = zeros(1,nnode);
%                 tempP = zeros(1,nline);
%                 tempQ = zeros(1,nline);
%                 tempu = zeros(1,nnode); tempu(:,k1) = 1;
%                 tempv = zeros(1,nnode);
%                 Aeq = [Aeq;
%                     zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
%                 beq = [beq;
%                     0];
%                 tempu = zeros(1,nnode);
%                 tempv = zeros(1,nnode); tempv(:,k1) = 1;
%                 Aeq = [Aeq;
%                     zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
%                 beq = [beq;
%                     0];
%             elseif wmaxpu(ph,k1) > 0
%                 tempE = zeros(ncirc+1,nnode);
%                 tempT = zeros(ncirc+1,nnode);
%                 tempP = zeros(ncirc+1,nline);
%                 tempQ = zeros(ncirc+1,nline);
%                 tempu = zeros(ncirc+1,nnode); tempu(:,k1) = cos(delta)';
%                 tempv = zeros(ncirc+1,nnode); tempv(:,k1) = sin(delta)';
%                 Aineq = [Aineq;
%                     zeros(ncirc+1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(ncirc+1,(3-ph)*nvar)];
%                 bineq = [bineq;
%                     wmaxpu(ph,k1)*ones(ncirc+1,1)];
%             end
%         end
%     end
% end
% 
% clear tempE tempT tempP tempQ tempu tempu tempv

end