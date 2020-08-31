function [nvar, Aineq, bineq, Aeq, beq] = create_CVX_Matrix_MagnitudeBalance_20180101(network, sim)
%%

% Michael Sankur - msankur@lbl.gov
% 2018.01.01

% This function creates the matrices and vectors for the inequality and
% equality constrainsts of the 

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
% sim - struct containing parameters pertinent to the simulation, such as
% objective function coefficients, higher order terms obtained by solving
% power flow

% OUTPUT(S)
% nvar - number of variables per phase
% Aineq - matrix for all inequality constraints Aineq X <= bineq
% bineq - RHS vector for all inequality constrains Aineq X <= bineq
% Aeq - matrix for all equality constraints Aeq X = beq
% beq - RHS vector for all equality constrains Aeq X = beq

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

% line paramters
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
rho = sim.rho;

% iteration parameters
if isempty(sim.VNR)
    VNR = [1*ones(1,nnode);
        1*exp(j*-120*pi/180)*ones(1,nnode);
        1*exp(j*120*pi/180)*ones(1,nnode)].*nodes.PH;
else
    VNR = sim.VNR;
end
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


%% Set up optimization matrices

Aineq = []; % Inequality constraint matrix
bineq = []; % Inequality constraint vector
Aeq = []; % Equality constraint matrix
beq = []; % Equality constraint vector

%% Voltage magnitude equality and inequality constraints

% If no phase present at node - Set voltage magnitude to zero
% If phase present at node - Add inequality constraint on squared voltage
% magnitude - 0.95^2 <= E_n^phi <= 1.05^2

for ph = 1:3
    for k1 = 2:nnode
        if nodes.PH(ph,k1) == 0
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
        elseif nodes.PH(ph,k1) == 1
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

% h = 1000
% ypen = 1/h*exp(-h*(x-0.95)) + 1/h*exp(h*(x-1.05));

%% Feeder head voltage magnitude (initialize) constraints

% Initialize feeder head square voltage magnitude as 1 for all three phases
% E_1 = [1 1 1]

tempE = zeros(1,nnode); tempE(1,1) = 1;
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
    abs(sim.Vnom)];

clear tempE tempT tempP tempQ tempu tempv

%% Feeder head voltage angle (initialize) constraints

% Initialize feeder head square voltage magnitude as 1 for all three phases
% theta_1^a = 0; theta_1^b = -120; theta_1^c = 120

tempE = zeros(1,nnode);
tempT = zeros(1,nnode); tempT(1,1) = 1;
tempP = zeros(1,nline);
tempQ = zeros(1,nline);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

Aeq = [Aeq;
    tempE tempT tempP tempQ tempu tempv, zeros(1,nvar), zeros(1,nvar);
    zeros(1,nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,nvar);
    zeros(1,nvar), zeros(1,nvar), tempE tempT tempP tempQ tempu tempv];

beq = [beq;
    180/pi*angle(sim.Vnom)];

clear tempE tempT tempP tempQ tempu tempv


%% Approximate |V| equality constraints

% |V_n^phi| = sqrt(y_n^phi)
% Approximate |V_n^phi| using first order Taylor Expansion
% |V_n^phi| ~ 0.5*(1 + y_n^phi)
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
%             tempVmag = zeros(1,n); tempVmag(k1) = 1;
%             
%             Aeq = [Aeq;
%                 zeros(1,(ph-1)*nvar), tempE tempVmag zeros(1,5*n), zeros(1,(3-ph)*nvar)];
%             beq = [beq;
%                 0.5];
%         end
%     end
% end

%% Power flow equality constraints

% If no phase present at node - Set power entering nodes to zero -
% P_n^phi = 0, Q_n^phi = 0
% If phase present at node -
% P_m^phi - u_m^phi - sum_n ( P_n^phi ) = p_n^phi + sum_n Re{L_mn^phi}
% Q_m^phi - v_m^phi - sum_n ( Q_n^phi ) = q_n^phi - cap_n^phi + sum_n Im{L_mn^phi}
% p_n^phi = p_n^phi*(A_PQ,n^phi + A_I,n^phi |V_n^phi| + A_Z,n^phi E_n^phi) + u_n^phi
% q_n^phi = q_n^phi*(A_PQ,n^phi + A_I,n^phi |V_n^phi| + A_Z,n^phi E_n^phi) - cap_n^phi + v_n^phi

tempE = zeros(1,nnode);
tempT = zeros(1,nnode);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

for ph = 1:3
    for k1 = 1:nline
        if lines.PH(ph,k1) == 0
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
        if nodes.PH(ph,k1) == 0
            for k2 = 1:size(nodes.inmat,1)
                if nodes.inmat(k2,k1) ~= 0
                    tempE = zeros(1,nnode);
                    tempT = zeros(1,nnode);
                    tempPQ = zeros(1,nline); tempPQ(1,nodes.inmat(k2,k1)) = 1;
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
            for k2 = 1:size(nodes.outmat,1)
                if nodes.outmat(k2,k1) ~= 0
                    tempE = zeros(1,nnode);
                    tempT = zeros(1,nnode);
                    tempPQ = zeros(1,nline); tempPQ(1,nodes.outmat(k2,k1)) = 1;
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
        elseif nodes.PH(ph,k1) == 1
            tempE = zeros(1,nnode); tempE(1,k1) = -spu(ph,k1)*aZ(ph,k1);
            % tempVmag = zeros(1,nnode); tempVmag(1,k1) = -spu(ph,k1)*aI(ph,k1);
            tempT = zeros(1,nnode);
            tempPQ = zeros(1,nline);
            tempPQ(nodes.inmat(nodes.inmat(:,k1) ~= 0,k1)) = 1;
            tempPQ(nodes.outmat(nodes.outmat(:,k1) ~= 0,k1)) = -1;
            tempuv = zeros(1,nnode); tempuv(1,k1) = -1;
            Aeq  = [Aeq;
                zeros(1,(ph-1)*nvar), real(tempE) tempT tempPQ zeros(1,nline) tempuv zeros(1,nnode), zeros(1,(3-ph)*nvar);
                zeros(1,(ph-1)*nvar), imag(tempE) tempT zeros(1,1*nline) tempPQ zeros(1,nnode) tempuv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                real(spu(ph,k1)*aPQ(ph,k1));
                imag(spu(ph,k1)*aPQ(ph,k1)) - cappu(ph,k1)];
%             for k2 = 1:nnode
%                 if feeder.FM(k1,k2) == 1
%                     beq(end-1:end) = beq(end-1:end) + [real(Lmn(ph,k2)); imag(Lmn(ph,k2))];
%                 end
%             end
        end        
    end        
end

clear tempE tempT tempP tempQ tempPQ tempu tempv tempuv

%% Calculate gamma terms for each node

for k1 = 1:nline
    txnode = TXnum(k1);
    gamman = VNR(:,txnode)*(1./VNR(:,txnode)).';
    if nodes.PH(1,txnode) == 0
        gamman(1,:) = 0; gamman(:,1) = 0;
    end
    if nodes.PH(2,txnode) == 0
        gamman(2,:) = 0; gamman(:,2) = 0;
    end
    if nodes.PH(3,txnode) == 0
        gamman(3,:) = 0; gamman(:,3) = 0;
    end
    gammaFZpu(:,:,k1) = gamman.*conj(FZpu(:,:,k1));
end
Mmn = real(gammaFZpu);
Nmn = -imag(gammaFZpu);

%% Voltage magnitude approximation equality constraints

% If no phase present at node - Set voltage at phase to zero - y_n^phi = 0
% If phase present at node - Propagate voltage up feeder (toward head)
% -E_m + E_n + M_mn P_n + N_mn Q_n = -H_mn

tempT = zeros(1,nnode);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

% Phase A
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if nodes.PH(1,txnode) == 1 && nodes.PH(1,rxnode) == 1 && lines.PH(1,k1) == 1
        tempE = zeros(1,nnode); tempE(1,txnode) = -1; tempE(1,rxnode) = 1;
        tempPa = zeros(1,nline); tempPa(1,k1) = 2*Mmn(1,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = 2*Mmn(1,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = 2*Mmn(1,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = 2*Nmn(1,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = 2*Nmn(1,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = 2*Nmn(1,3,k1);

        Aeq = [Aeq;
            tempE tempT tempPa tempQa tempu tempv, ...
            zeros(1,nnode) tempT tempPb tempQb tempu tempv, ...
            zeros(1,nnode) tempT tempPc tempQc tempu tempv];
        beq = [beq;
            0];
%             Hmn(1,k1)];
    end      
end

% Phase B
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);
    
    if nodes.PH(2,txnode) == 1 && nodes.PH(2,rxnode) == 1 && lines.PH(2,k1) == 1
        tempE = zeros(1,nnode); tempE(1,txnode) = -1; tempE(1,rxnode) = 1;
        tempPa = zeros(1,nline); tempPa(1,k1) = 2*Mmn(2,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = 2*Mmn(2,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = 2*Mmn(2,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = 2*Nmn(2,1,k1);   
        tempQb = zeros(1,nline); tempQb(1,k1) = 2*Nmn(2,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = 2*Nmn(2,3,k1);

        Aeq = [Aeq;
            zeros(1,nnode) tempT tempPa tempQa tempu tempv, ...
            tempE tempT tempPb tempQb tempu tempv, ...
            zeros(1,nnode) tempT tempPc tempQc tempu tempv];
        beq = [beq;
            0];
%             Hmn(1,k1)];
    end      
end

% Phase C
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if nodes.PH(3,txnode) == 1 && nodes.PH(3,rxnode) == 1 && lines.PH(3,k1) == 1
        tempE = zeros(1,nnode); tempE(1,txnode) = -1; tempE(1,rxnode) = 1;
        tempPa = zeros(1,nline); tempPa(1,k1) = 2*Mmn(3,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = 2*Mmn(3,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = 2*Mmn(3,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = 2*Nmn(3,1,k1);   
        tempQb = zeros(1,nline); tempQb(1,k1) = 2*Nmn(3,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = 2*Nmn(3,3,k1);

        Aeq = [Aeq;
            zeros(1,nnode) tempT tempPa tempQa tempu tempv, ...
            zeros(1,nnode) tempT tempPb tempQb tempu tempv, ...
            tempE tempT tempPc tempQc tempu tempv];
        beq = [beq;
            0];
%             Hmn(1,k1)];
    end      
end

clear tempE tempT tempu tempv
clear tempPa tempQa tempPb tempQb tempPc tempQc

%% Voltage angle approximation equality constraints

% If no phase present at node - No constraint on phase angle
% If phase present at node - Propagate voltage up feeder (toward head)
% |V_m||V_n|(-Theta_m + Theta_n) + 1/2 N_mn P_n - 1/2 M_mn Q_n = 0
% Above is not exact equation, read papers to find it

tempE = zeros(1,nnode);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

% Phase A
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if nodes.PH(1,txnode) == 1 && nodes.PH(1,rxnode) == 1 && lines.PH(1,k1) == 1
        tempT = zeros(1,nnode);
        tempT(1,txnode) = -1; tempT(1,rxnode) = 1;
        tempT = tempT*pi/180;
        tempPa = zeros(1,nline); tempPa(1,k1) = Nmn(1,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = Nmn(1,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = Nmn(1,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = -Mmn(1,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = -Mmn(1,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = -Mmn(1,3,k1);

        Aeq = [Aeq;
            tempE tempT tempPa tempQa tempu tempv, ...
            tempE zeros(1,nnode) tempPb tempQb tempu tempv, ...
            tempE zeros(1,nnode) tempPc tempQc tempu tempv];
        beq = [beq;
            0];
    end      
end

% Phase B
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if nodes.PH(2,txnode) == 1 && nodes.PH(2,rxnode) == 1 && lines.PH(2,k1) == 1
        tempT = zeros(1,nnode);
        tempT(1,txnode) = -1; tempT(1,rxnode) = 1;
        tempT = tempT*pi/180;
        tempPa = zeros(1,nline); tempPa(1,k1) = Nmn(2,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = Nmn(2,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = Nmn(2,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = -Mmn(2,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = -Mmn(2,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = -Mmn(2,3,k1);

        Aeq = [Aeq;
            tempE zeros(1,nnode) tempPa tempQa tempu tempv, ...
            tempE tempT tempPb tempQb tempu tempv, ...
            tempE zeros(1,nnode) tempPc tempQc tempu tempv];
        beq = [beq;
            0];
    end      
end

% Phase C
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if nodes.PH(3,txnode) == 1 && nodes.PH(3,rxnode) == 1 && lines.PH(3,k1) == 1
        tempT = zeros(1,nnode);
        tempT(1,txnode) = -1; tempT(1,rxnode) = 1;
        tempT = tempT*pi/180;
        tempPa = zeros(1,nline); tempPa(1,k1) = Nmn(3,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = Nmn(3,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = Nmn(3,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = -Mmn(3,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = -Mmn(3,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = -Mmn(3,3,k1);

        Aeq = [Aeq;
            tempE zeros(1,nnode) tempPa tempQa tempu tempv, ...
            tempE zeros(1,nnode) tempPb tempQb tempu tempv, ...
            tempE tempT tempPc tempQc tempu tempv];
        beq = [beq;
            0];
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
        if nodes.PH(ph,k1) == 0
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

if isfield(cons,'negativerealpower') == 1
    for ph = 1:3
        for k1 = 1:nnode
            if cons.negativerealpower(ph,k1) == 1
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
% cos(x) u_n^phi + sin(x) v_n^phi <= wmaxpu_n^phi
% x = 0:dx:2*pi

ncirc = 3600;
delta = linspace(0,2*pi,ncirc+1);

if isfield(cons,'wmaxpu') == 1
    for ph = 1:3
        for k1 = 1:nnode
            if wmaxpu(ph,k1) == 0
                tempE = zeros(1,nnode);
                tempT = zeros(1,nnode);
                tempP = zeros(1,nline);
                tempQ = zeros(1,nline);
                tempu = zeros(1,nnode); tempu(:,k1) = 1;
                tempv = zeros(1,nnode);
                Aeq = [Aeq;
                    zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
                beq = [beq;
                    0];
                tempE = zeros(1,nnode);
                tempT = zeros(1,nnode);
                tempP = zeros(1,nline);
                tempQ = zeros(1,nline);
                tempu = zeros(1,nnode);
                tempv = zeros(1,nnode); tempv(:,k1) = 1;
                Aeq = [Aeq;
                    zeros(1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
                beq = [beq;
                    0];
            elseif wmaxpu(ph,k1) > 0
                tempE = zeros(ncirc+1,nnode);
                tempT = zeros(ncirc+1,nnode);
                tempP = zeros(ncirc+1,nline);
                tempQ = zeros(ncirc+1,nline);
                tempu = zeros(ncirc+1,nnode); tempu(:,k1) = cos(delta)';
                tempv = zeros(ncirc+1,nnode); tempv(:,k1) = sin(delta)';
                Aineq = [Aineq;
                    zeros(ncirc+1,(ph-1)*nvar), tempE tempT tempP tempQ tempu tempv, zeros(ncirc+1,(3-ph)*nvar)];
                bineq = [bineq;
                    wmaxpu(ph,k1)*ones(ncirc+1,1)];
            end
        end
    end
end

clear tempE tempT tempP tempQ tempu tempu tempv

end