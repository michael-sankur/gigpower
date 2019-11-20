function OPTSOL = Solver_LinDist3Flow_TX_minP0_20160603_mag(feeder, nodes, lines, configs, loads, caps, controllers, sim)
%%

% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

% Variables : [Y_a; P_a; Q_a; u_a; v_a;
%               Y_b; P_b; Q_b; u_b; v_b;
%               Y_c; P_c; Q_c; u_c; v_c];

% Y_a = [y_1^a y_2^a ... y_N^a]

%%

% feeder parameters
% FM = feeder.FM;
% PH = feeder.PH;
% FZpu = feeder.FZpu;

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

% capacitor paramters

cappu = caps.cappu;

% controller parameters
wmaxpu = controllers.wmaxpu;

% sim parameters

rho = sim.rho

% iteration parameters
if isempty(sim.Vfbs)
    Vfbs = [1*ones(1,nnode);
        1*exp(j*-120*pi/180)*ones(1,nnode);
        1*exp(j*120*pi/180)*ones(1,nnode)].*nodes.PH;
else
    Vfbs = sim.Vfbs;
end
if isempty(sim.Lfbs)
    Lfbs = zeros(3,nline);
else
    Lfbs = sim.Lfbs;
end
if isempty(sim.Hfbs)
    Hfbs = zeros(3,nline);
else
    Hfbs = sim.Hfbs;
end

Vmagfbs = abs(Vfbs);

%%

nnode = nodes.nnode;
nline = lines.nline;

nvar = nnode + nline + nline + nnode + nnode; % number of variables per phase

%% Set up optimization matrices

Aineq = []; % Inequality constraint matrix
bineq = []; % Inequality constraint vector
Aeq = []; % Equality constraint matrix
beq = []; % Equality constraint vector

%% Voltage magnitude equality and inequality constraints

% If no phase present at node - Set voltage magnitude to zero
% If phase present at node - Add inequality constraint on squared voltage
% magnitude - 0.95^2 <= y_n^phi <= 1.05^2

for ph = 1:3
    for k1 = 2:nnode
        if nodes.PH(ph,k1) == 0
            tempY = zeros(1,nnode); tempY(1,k1) = 1;
            tempP = zeros(1,nline);
            tempQ = zeros(1,nline);
            tempu = zeros(1,nnode);
            tempv = zeros(1,nnode);
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), tempY tempP tempQ tempu tempv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0];
        elseif nodes.PH(ph,k1) == 1
            tempY = zeros(2,nnode); tempY(:,k1) = [1; -1];
            tempP = zeros(2,nline);
            tempQ = zeros(2,nline);
            tempu = zeros(2,nnode);
            tempv = zeros(2,nnode);
            Aineq = [Aineq;
                zeros(2,(ph-1)*nvar), tempY tempP tempQ tempu tempv, zeros(2,(3-ph)*nvar)];
            bineq = [bineq;
                1.05^2;
                -(0.95^2)];
        end
    end
end

clear tempY tempDelta tempP tempQ tempu tempv

% h = 1000
% ypen = 1/h*exp(-h*(x-0.95)) + 1/h*exp(h*(x-1.05));

%% Feeder head voltage magnitude (initialize) constraints

% Initialize feeder head square voltage magnitude as 1 for all three phases
% Y_1 = [1 1 1]

tempY = zeros(1,nnode); tempY(1,2) = 1;
tempDelta = zeros(1,nnode);
tempP = zeros(1,nline);
tempQ = zeros(1,nline);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

Aeq = [Aeq;
    tempY tempP tempQ tempu tempv, zeros(1,nvar), zeros(1,nvar);
    zeros(1,nvar), tempY tempP tempQ tempu tempv, zeros(1,nvar);
    zeros(1,nvar), zeros(1,nvar), tempY tempP tempQ tempu tempv];

beq = [beq;
    1;
    1;
    1];

clear tempY tempDelta tempP tempQ tempu tempv

%% Feeder head voltage angle (initialize) constraints

% Initialize feeder head square voltage magnitude as 1 for all three phases
% theta_1^a = 0; theta_1^b = -120; theta_1^c = 120

% tempY = zeros(1,nnode);
% tempDelta = zeros(1,nnode); tempDelta(1,2) = 1;
% tempP = zeros(1,nline);
% tempQ = zeros(1,nline);
% tempu = zeros(1,nnode);
% tempv = zeros(1,nnode);
% 
% Aeq = [Aeq;
%     tempY tempP tempQ tempu tempv, zeros(1,nvar), zeros(1,nvar);
%     zeros(1,nvar), tempY tempP tempQ tempu tempv, zeros(1,nvar);
%     zeros(1,nvar), zeros(1,nvar), tempY tempP tempQ tempu tempv];
% 
% beq = [beq;
%     0;
%     -120;
%     120];
% 
% clear tempY tempDelta tempP tempQ tempu tempv


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
%             tempY = zeros(1,n); tempY(k1) = -0.5;
%             tempVmag = zeros(1,n); tempVmag(k1) = 1;
%             
%             Aeq = [Aeq;
%                 zeros(1,(ph-1)*nvar), tempY tempVmag zeros(1,5*n), zeros(1,(3-ph)*nvar)];
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
% p_n^phi = p_n^phi*(A_PQ,n^phi + A_I,n^phi |V_n^phi| + A_Z,n^phi y_n^phi) + u_n^phi
% q_n^phi = q_n^phi*(A_PQ,n^phi + A_I,n^phi |V_n^phi| + A_Z,n^phi y_n^phi) -cap_n^phi + v_n^phi

tempY = zeros(1,nnode);
tempDelta = zeros(1,nnode);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

for ph = 1:3
    for k1 = 1:nline
        if lines.PH(ph,k1) == 0
            tempPQ = zeros(1,nline); tempPQ(1,k1) = 1;
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), tempY tempPQ zeros(1,nline) tempu tempv, zeros(1,(3-ph)*nvar);
                zeros(1,(ph-1)*nvar), tempY zeros(1,nline) tempPQ tempu tempv, zeros(1,(3-ph)*nvar)];
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
                    tempY = zeros(1,nnode);
                    tempDelta = zeros(1,nnode);
                    tempPQ = zeros(1,nline); tempPQ(1,nodes.inmat(k2,k1)) = 1;
                    tempu = zeros(1,nnode);
                    tempv = zeros(1,nnode);                    
                    Aeq = [Aeq;
                        zeros(1,(ph-1)*nvar), tempY tempPQ zeros(1,nline) tempu tempv, zeros(1,(3-ph)*nvar);
                        zeros(1,(ph-1)*nvar), tempY zeros(1,nline) tempPQ tempu tempv, zeros(1,(3-ph)*nvar)];
                    beq = [beq;
                        0;
                        0];                   
                end
            end
            for k2 = 1:size(nodes.outmat,1)
                if nodes.outmat(k2,k1) ~= 0
                    tempY = zeros(1,nnode);
                    tempDelta = zeros(1,nnode);
                    tempPQ = zeros(1,nline); tempPQ(1,nodes.outmat(k2,k1)) = 1;
                    tempu = zeros(1,nnode);
                    tempv = zeros(1,nnode);
                    Aeq = [Aeq;
                            zeros(1,(ph-1)*nvar), tempY tempPQ zeros(1,nline) tempu tempv, zeros(1,(3-ph)*nvar);
                            zeros(1,(ph-1)*nvar), tempY zeros(1,nline) tempPQ tempu tempv, zeros(1,(3-ph)*nvar)];
                    beq = [beq;
                        0;
                        0];
                end
            end
        elseif nodes.PH(ph,k1) == 1
            tempY = zeros(1,nnode); tempY(1,k1) = -spu(ph,k1)*aZ(ph,k1);
            % tempVmag = zeros(1,nnode); tempVmag(1,k1) = -spu(ph,k1)*aI(ph,k1);
            tempDelta = zeros(1,nnode);
            tempPQ = zeros(1,nline);
            tempPQ(nodes.inmat(nodes.inmat(:,k1) ~= 0,k1)) = 1;
            tempPQ(nodes.outmat(nodes.outmat(:,k1) ~= 0,k1)) = -1;
            tempuv = zeros(1,nnode); tempuv(1,k1) = -1;
            Aeq  = [Aeq;
                zeros(1,(ph-1)*nvar), real(tempY) tempPQ zeros(1,nline) tempuv zeros(1,nnode), zeros(1,(3-ph)*nvar);
                zeros(1,(ph-1)*nvar), imag(tempY) zeros(1,1*nline) tempPQ zeros(1,nnode) tempuv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                real(spu(ph,k1)*aPQ(ph,k1));
                imag(spu(ph,k1)*aPQ(ph,k1)) - cappu(ph,k1)];
%             for k2 = 1:nnode
%                 if feeder.FM(k1,k2) == 1
%                     beq(end-1:end) = beq(end-1:end) + [real(Lfbs(ph,k2)); imag(Lfbs(ph,k2))];
%                 end
%             end
        end        
    end        
end

clear tempY tempDelta tempP tempQ tempPQ tempu tempv tempuv

%% Calculate gamma terms for each node

for k1 = 1:nline
    txnode = TXnum(k1);
    gamman = Vfbs(:,txnode)*(1./Vfbs(:,txnode)).';
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
% -Y_m + Y_n + M_mn P_n + N_mn Q_n = -H_mn

tempDelta = zeros(1,nnode);
tempu = zeros(1,nnode);
tempv = zeros(1,nnode);

% Phase A
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if nodes.PH(1,txnode) == 1 && nodes.PH(1,rxnode) == 1 && lines.PH(1,k1) == 1
        tempY = zeros(1,nnode); tempY(1,txnode) = -1; tempY(1,rxnode) = 1;
        tempPa = zeros(1,nline); tempPa(1,k1) = 2*Mmn(1,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = 2*Mmn(1,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = 2*Mmn(1,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = 2*Nmn(1,1,k1);
        tempQb = zeros(1,nline); tempQb(1,k1) = 2*Nmn(1,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = 2*Nmn(1,3,k1);

        Aeq = [Aeq;
            tempY tempPa tempQa tempu tempv, ...
            zeros(1,nnode) tempPb tempQb tempu tempv, ...
            zeros(1,nnode) tempPc tempQc tempu tempv];
        beq = [beq;
            0];
%             Hfbs(1,k1)];
    end      
end

% Phase B
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);
    
    if nodes.PH(2,txnode) == 1 && nodes.PH(2,rxnode) == 1 && lines.PH(2,k1) == 1
        tempY = zeros(1,nnode); tempY(1,txnode) = -1; tempY(1,rxnode) = 1;
        tempPa = zeros(1,nline); tempPa(1,k1) = 2*Mmn(2,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = 2*Mmn(2,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = 2*Mmn(2,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = 2*Nmn(2,1,k1);   
        tempQb = zeros(1,nline); tempQb(1,k1) = 2*Nmn(2,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = 2*Nmn(2,3,k1);

        Aeq = [Aeq;
            zeros(1,nnode) tempPa tempQa tempu tempv, ...
            tempY tempPb tempQb tempu tempv, ...
            zeros(1,nnode) tempPc tempQc tempu tempv];
        beq = [beq;
            0];
%             Hfbs(1,k1)];
    end      
end

% Phase C
for k1 = 1:nline
    
    txnode = TXnum(k1);
    rxnode = RXnum(k1);

    if nodes.PH(3,txnode) == 1 && nodes.PH(3,rxnode) == 1 && lines.PH(3,k1) == 1
        tempY = zeros(1,nnode); tempY(1,txnode) = -1; tempY(1,rxnode) = 1;
        tempPa = zeros(1,nline); tempPa(1,k1) = 2*Mmn(3,1,k1);
        tempPb = zeros(1,nline); tempPb(1,k1) = 2*Mmn(3,2,k1);
        tempPc = zeros(1,nline); tempPc(1,k1) = 2*Mmn(3,3,k1);
        tempQa = zeros(1,nline); tempQa(1,k1) = 2*Nmn(3,1,k1);   
        tempQb = zeros(1,nline); tempQb(1,k1) = 2*Nmn(3,2,k1);
        tempQc = zeros(1,nline); tempQc(1,k1) = 2*Nmn(3,3,k1);

        Aeq = [Aeq;
            zeros(1,nnode) tempPa tempQa tempu tempv, ...
            zeros(1,nnode) tempPb tempQb tempu tempv, ...
            tempY tempPc tempQc tempu tempv];
        beq = [beq;
            0];
%             Hfbs(1,k1)];
    end      
end

clear tempY tempDelta tempu tempv
clear tempPa tempQa tempPb tempQb tempPc tempQc

%% Voltage angle approximation equality constraints

% If no phase present at node - No constraint on phase angle
% If phase present at node - Propagate voltage up feeder (toward head)
% |V_m||V_n|(-Theta_m + Theta_n) + 1/2 N_mn P_n - 1/2 M_mn Q_n = 0
% Above is not exact equation, read papers to find it

% tempY = zeros(1,nnode);
% tempu = zeros(1,nnode);
% tempv = zeros(1,nnode);
% 
% % Phase A
% for k1 = 1:nline
%     
%     txnode = TXnum(k1);
%     rxnode = RXnum(k1);
% 
%     if nodes.PH(1,txnode) == 1 && nodes.PH(1,rxnode) == 1 && lines.PH(1,k1) == 1
%         tempDelta = zeros(1,nnode);
%         tempDelta(1,txnode) = -1; tempDelta(1,rxnode) = 1;
%         tempDelta = tempDelta*pi/180;
%         tempPa = zeros(1,nline); tempPa(1,k1) = Nmn(1,1,k1);
%         tempPb = zeros(1,nline); tempPb(1,k1) = Nmn(1,2,k1);
%         tempPc = zeros(1,nline); tempPc(1,k1) = Nmn(1,3,k1);
%         tempQa = zeros(1,nline); tempQa(1,k1) = -Mmn(1,1,k1);
%         tempQb = zeros(1,nline); tempQb(1,k1) = -Mmn(1,2,k1);
%         tempQc = zeros(1,nline); tempQc(1,k1) = -Mmn(1,3,k1);
% 
%         Aeq = [Aeq;
%             tempY tempDelta tempPa tempQa tempu tempv, ...
%             tempY zeros(1,nnode) tempPb tempQb tempu tempv, ...
%             tempY zeros(1,nnode) tempPc tempQc tempu tempv];
%         beq = [beq;
%             0];
%     end      
% end
% 
% % Phase B
% for k1 = 1:nline
%     
%     txnode = TXnum(k1);
%     rxnode = RXnum(k1);
% 
%     if nodes.PH(2,txnode) == 1 && nodes.PH(2,rxnode) == 1 && lines.PH(2,k1) == 1
%         tempDelta = zeros(1,nnode);
%         tempDelta(1,txnode) = -1; tempDelta(1,rxnode) = 1;
%         tempDelta = tempDelta*pi/180;
%         tempPa = zeros(1,nline); tempPa(1,k1) = Nmn(2,1,k1);
%         tempPb = zeros(1,nline); tempPb(1,k1) = Nmn(2,2,k1);
%         tempPc = zeros(1,nline); tempPc(1,k1) = Nmn(2,3,k1);
%         tempQa = zeros(1,nline); tempQa(1,k1) = -Mmn(2,1,k1);
%         tempQb = zeros(1,nline); tempQb(1,k1) = -Mmn(2,2,k1);
%         tempQc = zeros(1,nline); tempQc(1,k1) = -Mmn(2,3,k1);
% 
%         Aeq = [Aeq;
%             tempY zeros(1,nnode) tempPa tempQa tempu tempv, ...
%             tempY tempDelta tempPb tempQb tempu tempv, ...
%             tempY zeros(1,nnode) tempPc tempQc tempu tempv];
%         beq = [beq;
%             0];
%     end      
% end
% 
% % Phase C
% for k1 = 1:nline
%     
%     txnode = TXnum(k1);
%     rxnode = RXnum(k1);
% 
%     if nodes.PH(3,txnode) == 1 && nodes.PH(3,rxnode) == 1 && lines.PH(3,k1) == 1
%         tempDelta = zeros(1,nnode);
%         tempDelta(1,txnode) = -1; tempDelta(1,rxnode) = 1;
%         tempDelta = tempDelta*pi/180;
%         tempPa = zeros(1,nline); tempPa(1,k1) = Nmn(3,1,k1);
%         tempPb = zeros(1,nline); tempPb(1,k1) = Nmn(3,2,k1);
%         tempPc = zeros(1,nline); tempPc(1,k1) = Nmn(3,3,k1);
%         tempQa = zeros(1,nline); tempQa(1,k1) = -Mmn(3,1,k1);
%         tempQb = zeros(1,nline); tempQb(1,k1) = -Mmn(3,2,k1);
%         tempQc = zeros(1,nline); tempQc(1,k1) = -Mmn(3,3,k1);
% 
%         Aeq = [Aeq;
%             tempY zeros(1,nnode) tempPa tempQa tempu tempv, ...
%             tempY zeros(1,nnode) tempPb tempQb tempu tempv, ...
%             tempY tempDelta tempPc tempQc tempu tempv];
%         beq = [beq;
%             0];
%     end      
% end
% 
% clear tempDelta tempPa tempQa tempPb tempQb tempPc tempQc

%% Zero inverter outputs for non control nodes and nonexistent phases

% For all phases at all nodes that are non control, set inverter output to
% zero - u_n^phi = 0, v_n^phi = 0

tempY = zeros(1,nnode);
tempDelta = zeros(1,nnode);
tempP = zeros(1,nline);
tempQ = zeros(1,nline);

for ph = 1:3
    for k1 = 2:nnode
        if nodes.PH(ph,k1) == 0
            tempuv = zeros(1,nnode); tempuv(1,k1) = 1;
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), tempY tempP tempQ tempuv zeros(1,nnode), zeros(1,(3-ph)*nvar);
                zeros(1,(ph-1)*nvar), tempY tempP tempQ zeros(1,nnode) tempuv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0;
                0];
        end
    end
end

clear tempY tempDelta tempP tempQ tempu tempu tempv

%% Zero inverter real power output for all DER

% For all phases at all nodes that are non control, set inverter output to
% zero - u_n^phi = 0, v_n^phi = 0

% for ph = 1:3
%     for k1 = 1:n
%         if cn(ph,k1) == 1
%             tempuv = zeros(1,n); tempuv(1,k1) = 1;
%             Aeq = [Aeq;
%                 zeros(1,(ph-1)*nvar), zeros(1,5*n) tempuv zeros(1,n), zeros(1,(3-ph)*nvar)];
%             beq = [beq;
%                 0];
%         end
%     end
% end
% 
% clear tempuv

%% Inverter control bounds - square

% for ph = 1:3
%     for k1 = 1:n
%         if cn(k1) == 1
%             if PH(ph,k1) == 0
%                 tempuv = zeros(1,n); tempuv(:,k1) = 1;
%                 Aeq = [Aeq;
%                     zeros(1,(ph-1)*nvar), zeros(1,3*n) tempuv zeros(1,n), zeros(1,(3-ph)*nvar);
%                     zeros(1,(ph-1)*nvar), zeros(1,3*n) zeros(1,n) tempuv, zeros(1,(3-ph)*nvar)];
%                 beq = [beq;
%                     0;
%                     0];
%             elseif PH(ph,k1) == 1
%                 tempuv = zeros(2,n); tempuv(:,k1) = [1; -1];
%                 Aineq = [Aineq;
%                     zeros(2,(ph-1)*nvar), zeros(2,3*n) tempuv zeros(2,n), zeros(2,(3-ph)*nvar);
%                     zeros(2,(ph-1)*nvar), zeros(2,3*n) zeros(2,n) tempuv, zeros(2,(3-ph)*nvar)];
%                 bineq = [bineq;
%                     wmax(k1);
%                     wmax(k1);
%                     wmax(k1);
%                     wmax(k1)]; 
%             end                       
%         end
%     end
% end
% 
% clear tempuv

%% Inverter control bounds - circle

% If no phase present at node - Set magnitude of inverter output to zero -
% |w_n^phi| <= 0
% If phase present at node - Set magnitude of inverter output to maximum
% for particular phase and node - |w_n^phi| <= wmax_n^phi

% Acirc = zeros(2,3*nvar,3*n);
% bcirc = zeros(n);
% 
% for ph = 1:3
%     for k1 = 1:n
%         Acirc(1, (ph-1)*nvar + 5*n + k1,(ph-1)*n + k1) = 1;
%         Acirc(2, (ph-1)*nvar + 6*n + k1,(ph-1)*n + k1) = 1;
%         if cn(ph,k1) == 1 && PH(ph,k1) == 0
%             bcirc((ph-1)*n + k1) = 0;
%         elseif  cn(ph,k1) == 1 && PH(ph,k1) == 1
%             bcirc((ph-1)*n + k1) = wmaxpu(ph,k1);  
%         end
%     end
% end

%% Inverter ramp rate constraints

% for ph = 1:3
%     for k1 = 1:n
%         if cn(ph,k1) == 1
%             tempuv = zeros(2,n); tempuv(:,k1) = [1; -1];
%             Aineq = [Aineq;
%                 zeros(2,(ph-1)*nvar), zeros(2,4*n) tempuv zeros(2,n), zeros(2,(3-ph)*nvar);
%                 zeros(2,(ph-1)*nvar), zeros(2,4*n) zeros(2,n) tempuv, zeros(2,(3-ph)*nvar)];
%             bineq = [bineq;
%                 rrpu(ph,k1) + real(wk1(ph,k1));
%                 rrpu(ph,k1) - real(wk1(ph,k1));
%                 rrpu(ph,k1) + imag(wk1(ph,k1));
%                 rrpu(ph,k1) - imag(wk1(ph,k1))];            
%         end
%     end
% end

% Aramp = zeros(2,3*nvar,3*n);
% bramp = zeros(3*n,1);
% 
% for ph = 1:3
%     for k1 = 1:n
%         Aramp(1,(ph-1)*nvar + 5*n + k1,(ph-1)*n + k1) = 1;
%         Aramp(2,(ph-1)*nvar + 6*n + k1,(ph-1)*n + k1) = 1;
%         if cn(ph,k1) == 1 && PH(ph,k1) == 0
%             bramp((ph-1)*n + k1) = 0;
%         elseif cn(ph,k1) == 1 && PH(ph,k1) == 1
%             bramp((ph-1)*n + k1) = rrpu(ph,k1);
%         end
%     end
% end
% 
% clear tempuv

%% Objective function - corrected

% Matrices to create a scalar or vector of voltage magnitude differences
% for each node. If a node has one phase, the vector is null. If the node
% has two phases, the vector contains Y_m,k - Y_n,k where m and n are the
% existing phases. If the node has three phases, the vector is as follows:
% [Y_a,k - Y_b,k ; Y_b,k - Y_c,k ; Y_a,k - Y_c,k].
% FV = [];
% for k1 = 2:n
%     tempY = zeros(1,n); tempY(k1) = 1;
%     if sum(PH(:,k1)) == 1
%         FV(:,:,k1) = zeros(3,3*nvar);
%     elseif PH(1,k1) == 1 && PH(2,k1) == 1 && PH(3,k1) == 0
%         FV(:,:,k1) = [tempY zeros(1,6*n), -tempY zeros(1,6*n), zeros(1,nvar);
%             zeros(2,nvar), zeros(2,nvar), zeros(2,nvar)];
%     elseif PH(1,k1) == 0 && PH(2,k1) == 1 && PH(3,k1) == 1
%         FV(:,:,k1) = [zeros(1,nvar), tempY zeros(1,6*n), -tempY zeros(1,6*n);
%             zeros(2,nvar), zeros(2,nvar), zeros(2,nvar)];
%     elseif PH(1,k1) == 1 && PH(2,k1) == 0 && PH(3,k1) == 1
%         FV(:,:,k1) = [tempY zeros(1,6*n), zeros(1,nvar), -tempY zeros(1,6*n);
%             zeros(2,nvar), zeros(2,nvar), zeros(2,nvar)];   
%     elseif PH(1,k1) == 1 && PH(2,k1) == 1 && PH(3,k1) == 1
%         FV(:,:,k1) = [tempY zeros(1,6*n), -tempY zeros(1,6*n), zeros(1,nvar);
%             zeros(1,nvar), tempY zeros(1,6*n), -tempY zeros(1,6*n);
%             tempY zeros(1,6*n), zeros(1,nvar), -tempY zeros(1,6*n)];
%     end
% end
% FV;
% size(FV);

% Matrix weighting control effort
% Fu1 = blkdiag(zeros(5*n,5*n),diag(controllers.cn(1,:)),diag(controllers.cn(1,:)));
% Fu2 = blkdiag(zeros(5*n,5*n),diag(controllers.cn(2,:)),diag(controllers.cn(2,:)));
% Fu3 = blkdiag(zeros(5*n,5*n),diag(controllers.cn(3,:)),diag(controllers.cn(3,:)));
% Fu = blkdiag(Fu1,Fu2,Fu3);

% Matrix weigting control effort
% FU = [];
% for k1 = 2:n
%     tempuv = zeros(1,n); tempuv(k1) = 1;
%     tempphase = [zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) tempuv zeros(1,n);
%         zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) tempuv];
%     FU(:,:,k1) = [tempphase zeros(2,nvar) zeros(2,nvar);
%         zeros(2,nvar) tempphase zeros(2,nvar);
%         zeros(2,nvar) zeros(2,nvar) tempphase];
% end

% Vector for feeder head power in objective function
fP0 = zeros(nvar,1);
fP0(nnode+1) = 1;
fP0 = [fP0; fP0; fP0];

fYA3ABC = zeros(1,nnode);
fYA3ABC(1,3) = 1;
fYA3ABC = [fYA3ABC zeros(1,nline) zeros(1,nline) zeros(1,nnode) zeros(1,nnode)]';
fYA3ABC = [1*fYA3ABC; 1*fYA3ABC; 1*fYA3ABC];

kYa = 3;
kYb = 3+nvar;
kYc = 3+2*nvar;

cvx_begin quiet
    expressions Z1 Z2;
    variable X(3*nvar);
%     Z1 = fP0'*X;
%     Z1 = fYA3ABC'*X;
%     Z1 = (X(kYa) - X(kYb))^2 + (X(kYb) - X(kYc))^2 + (X(kYa) - X(kYc))^2;
    for ph = 1:3
        for k1 = 2:nline
            kP = (ph-1)*nvar + nnode + k1;
            kQ = (ph-1)*nvar + nnode + nline + k1;
            Z1 = Z1 + norm(X([kP kQ]),2);
        end
    end
    minimize(Z1)
    subject to;
    Aeq * X == beq;
    Aineq * X <= bineq;
    for ph = 1:3
        for k1 = 1:nnode
            ku = (ph-1)*nvar + nnode +  nline + nline + k1;
            kv = (ph-1)*nvar + nnode + nline + nline + nnode + k1;
            norm(X([ku kv]),2) <= wmaxpu(ph,k1)
        end
    end
%     for ph = 1:3
%         for k1 = 2:nline
%             kP = (ph-1)*nvar + nnode + k1;
%             kQ = (ph-1)*nvar + nnode + nline + k1;
%             norm(X([kP kQ]),2) <= 1.5
%         end
%     end
cvx_end

OPTSOL.X = X;
OPTSOL.nvar = nvar;
OPTSOL.cvx_optval = cvx_optval;
OPTSOL.cvx_status = cvx_status;
OPTSOL.Z1 = Z1;
OPTSOL.Z2 = Z2;

end
