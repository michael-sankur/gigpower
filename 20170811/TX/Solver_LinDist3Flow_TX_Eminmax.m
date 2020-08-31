function [OPTSOLMIN OPTSOLMAX] = Solver_LinDist3Flow_TX_Eminmax(feeder,loads,controllers,sim)
%%

% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

% Variables : [Y_a; Vmag_a; D_a; P_a; Q_a; u_a; v_a;
%               Y_b; Vmag_b; D_b; P_b; Q_b; u_b; v_b;
%               Y_c; Vmag_c; D_c; P_c; Q_c; u_c; v_c];

%%

% feeder parameters
n = feeder.n;
FM = feeder.FM;
PH = feeder.PH;
FZpu = feeder.FZpu;

% load parameters
spu = loads.spu;
APQ = loads.APQ;
AI = loads.AI;
AZ = loads.AZ;
cappu = loads.cappu;

% controller parameters
cn = controllers.cn;
cnstate = controllers.cnstate;
wmaxpu = controllers.wmaxpu;
rrpu = controllers.rrpu;
wk1 = controllers.wk1;

% sim parameters

tn = sim.tn;
tp = sim.tp;

% iteration parameters
if isempty(sim.Vfbs)
    Vfbs = [1*ones(1,n);
        1*exp(j*-120*pi/180)*ones(1,n);
        1*exp(j*120*pi/180)*ones(1,n)].*PH;
end
if isempty(sim.Lfbs)
    Lfbs = zeros(3,n);
end
if isempty(sim.Hfbs)
    Hfbs = zeros(3,n);
end

Vmagfbs = abs(Vfbs);

%%

nvar = 7*n; % number of variables per phase

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
    for k1 = 2:n
        if PH(ph,k1) == 0
%             tempY = zeros(1,n); tempY(1,k1) = 1;
%             Aeq = [Aeq;
%                 zeros(1,(ph-1)*nvar), tempY zeros(1,6*n), zeros(1,(3-ph)*nvar)];
%             beq = [beq;
%                 0];
        elseif PH(ph,k1) == 1
            tempY = zeros(2,n); tempY(:,k1) = [1; -1];
            Aineq = [Aineq;
                zeros(2,(ph-1)*nvar), tempY zeros(2,6*n), zeros(2,(3-ph)*nvar)];
            bineq = [bineq;
                1.05^2;
                -(0.95^2)];
        end
    end
end

clear tempY

% h = 1000
% ypen = 1/h*exp(-h*(x-0.95)) + 1/h*exp(h*(x-1.05));

%% Feeder head voltage magnitude (initialize) constraints

% Initialize feeder head square voltage magnitude as 1 for all three phases
% Y_1 = [1 1 1]
tempY = zeros(1,n); tempY(1) = 1;

Aeq = [Aeq;
    tempY zeros(1,6*n), zeros(1,nvar) zeros(1,nvar);
    zeros(1,nvar), tempY zeros(1,6*n), zeros(1,nvar);
    zeros(1,nvar), zeros(1,nvar), tempY zeros(1,6*n)];

beq = [beq;
    1;
    1;
    1];

clear tempY

%% Feeder head voltage angle (initialize) constraints

% Initialize feeder head square voltage magnitude as 1 for all three phases
% theta_1^a = 0; theta_1^b = -120; theta_1^c = 120

tempDelta = zeros(1,n); tempDelta(1,1) = 1;

Aeq = [Aeq;
    zeros(1,n) zeros(1,n) tempDelta zeros(1,4*n), zeros(1,nvar) zeros(1,nvar);
    zeros(1,nvar), zeros(1,n) zeros(1,n) tempDelta zeros(1,4*n), zeros(1,nvar);
    zeros(1,nvar), zeros(1,nvar), zeros(1,n) zeros(1,n) tempDelta zeros(1,4*n)];

beq = [beq;
    0;
    -120;
    120];

%% Approximate |V| equality constraints

% |V_n^phi| = sqrt(y_n^phi)
% Approximate |V_n^phi| using first order Taylor Expansion
% |V_n^phi| ~ 0.5*(1 + y_n^phi)
for ph = 1:3
    for k1 = 1:n
        if PH(ph,k1) == 0
            tempVmag = zeros(1,n); tempVmag(k1) = 1;
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), zeros(1,n) tempVmag zeros(1,5*n), zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0];
        elseif PH(ph,k1) == 1
            tempY = zeros(1,n); tempY(k1) = -0.5;
            tempVmag = zeros(1,n); tempVmag(k1) = 1;
            
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), tempY tempVmag zeros(1,5*n), zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0.5];
        end
    end
end

%% Calculate gamma terms for each node

for k1 = 2:n
    gamma(:,:,k1) = zeros(3,3);
    if sum(PH(:,k1) == [1 0 0]') == 3
        gamma(1,1,k1) = 1;
    elseif sum(PH(:,k1) == [0 1 0]') == 3
        gamma(2,2,k1) = 1;
    elseif sum(PH(:,k1) == [0 0 1]') == 3
        gamma(3,3,k1) = 1;
    elseif sum(PH(:,k1) == [1 1 0]') == 3
        gamma(:,:,k1) = [1 Vfbs(1,k1)/Vfbs(2,k1) 0;
            Vfbs(2,k1)/Vfbs(1,k1) 1 0;
            0 0 0];
    elseif sum(PH(:,k1) == [0 1 1]') == 3
        gamma(:,:,k1) = [0 0 0;
            0 1 Vfbs(2,k1)/Vfbs(3,k1);
            0 Vfbs(3,k1)/Vfbs(2,k1) 1];
    elseif sum(PH(:,k1) == [1 0 1]') == 3
        gamma(:,:,k1) = [1 0 Vfbs(1,k1)/Vfbs(3,k1);
            0 0 0;
            Vfbs(3,k1)/Vfbs(1,k1) 0 1];
    elseif sum(PH(:,k1) == [1 1 1]') == 3
        gamma(:,:,k1) = [1 Vfbs(1,k1)/Vfbs(2,k1) Vfbs(1,k1)/Vfbs(3,k1);
            Vfbs(2,k1)/Vfbs(1,k1) 1 Vfbs(2,k1)/Vfbs(3,k1);
            Vfbs(3,k1)/Vfbs(1,k1) Vfbs(3,k1)/Vfbs(2,k1) 1];
    end
    % multiply gamma_n^{phi psi} and {Z_{mn}^{phi psi}}^*
    gammaFZpu(:,:,k1) = gamma(:,:,k1).*conj(FZpu(:,:,k1));
end
gamma;

%% Voltage magnitude approximation equality constraints

% If no phase present at node - Set voltage at phase to zero - y_n^phi = 0
% If phase present at node - Propagate voltage up feeder (toward head)
% -Y_m + Y_n + M_mn P_n + N_mn Q_n = H_mn

% Phase A
for k1 = 2:n
    if PH(1,k1) == 0
        tempY = zeros(1,n); tempY(1,k1) = 1;
        Aeq = [Aeq;
            tempY zeros(1,6*n), zeros(1,nvar), zeros(1,nvar)];
        beq = [beq;
            0];
    elseif PH(1,k1) == 1        
        tempY = -[FM(k1,:)==-1]; tempY(1,k1) = 1;        
        tempPa = zeros(1,n); tempPa(1,k1) = 2*real(gammaFZpu(1,1,k1));
        tempPb = zeros(1,n); tempPb(1,k1) = 2*real(gammaFZpu(1,2,k1));
        tempPc = zeros(1,n); tempPc(1,k1) = 2*real(gammaFZpu(1,3,k1));
        tempQa = zeros(1,n); tempQa(1,k1) = -2*imag(gammaFZpu(1,1,k1));   
        tempQb = zeros(1,n); tempQb(1,k1) = -2*imag(gammaFZpu(1,2,k1));
        tempQc = zeros(1,n); tempQc(1,k1) = -2*imag(gammaFZpu(1,3,k1));

        Aeq = [Aeq;
            tempY zeros(1,n) zeros(1,n) tempPa tempQa zeros(1,2*n), ...
            zeros(1,n) zeros(1,n) zeros(1,n) tempPb tempQb zeros(1,2*n), ...
            zeros(1,n) zeros(1,n) zeros(1,n) tempPc tempQc zeros(1,2*n)];
        beq = [beq;
            -Hfbs(1,k1)];
    end      
end

% Phase B
for k1 = 2:n
    if PH(2,k1) == 0
        tempY = zeros(1,n); tempY(1,k1) = 1;
        Aeq = [Aeq;
            zeros(1,nvar), tempY zeros(1,6*n), zeros(1,nvar)];
        beq = [beq;
            0];
    elseif PH(2,k1) == 1
        tempY = -[FM(k1,:)==-1]; tempY(1,k1) = 1;       
        tempPa = zeros(1,n); tempPa(1,k1) = 2*real(gammaFZpu(2,1,k1));
        tempPb = zeros(1,n); tempPb(1,k1) = 2*real(gammaFZpu(2,2,k1));
        tempPc = zeros(1,n); tempPc(1,k1) = 2*real(gammaFZpu(2,3,k1));
        tempQa = zeros(1,n); tempQa(1,k1) = -2*imag(gammaFZpu(2,1,k1));
        tempQb = zeros(1,n); tempQb(1,k1) = -2*imag(gammaFZpu(2,2,k1));
        tempQc = zeros(1,n); tempQc(1,k1) = -2*imag(gammaFZpu(2,3,k1));

        Aeq = [Aeq;
            zeros(1,n) zeros(1,n) zeros(1,n) tempPa tempQa zeros(1,2*n), ...
            tempY zeros(1,n) zeros(1,n) tempPb tempQb zeros(1,2*n), ...
            zeros(1,n) zeros(1,n) zeros(1,n) tempPc tempQc zeros(1,2*n)];
        beq = [beq;
            -Hfbs(2,k1)];
    end
end

% Phase C
for k1 = 2:n
    if PH(3,k1) ==0
        tempY = zeros(1,n); tempY(1,k1) = 1;
        Aeq = [Aeq;
            zeros(1,nvar), zeros(1,nvar), tempY zeros(1,6*n)];
        beq = [beq;
            0];
    elseif PH(3,k1) == 1
        tempY = -[FM(k1,:)==-1]; tempY(1,k1) = 1;        
        tempPa = zeros(1,n); tempPa(1,k1) = 2*real(gammaFZpu(3,1,k1));
        tempPb = zeros(1,n); tempPb(1,k1) = 2*real(gammaFZpu(3,2,k1));
        tempPc = zeros(1,n); tempPc(1,k1) = 2*real(gammaFZpu(3,3,k1));
        tempQa = zeros(1,n); tempQa(1,k1) = -2*imag(gammaFZpu(3,1,k1));
        tempQb = zeros(1,n); tempQb(1,k1) = -2*imag(gammaFZpu(3,2,k1));
        tempQc = zeros(1,n); tempQc(1,k1) = -2*imag(gammaFZpu(3,3,k1));

        Aeq = [Aeq;
            zeros(1,n) zeros(1,n) zeros(1,n) tempPa tempQa zeros(1,2*n), ...
            zeros(1,n) zeros(1,n) zeros(1,n) tempPb tempQb zeros(1,2*n), ...
            tempY zeros(1,n) zeros(1,n) tempPc tempQc zeros(1,2*n)];
        beq = [beq;
            -Hfbs(3,k1)];
    end
end

clear tempY tempPa tempQa tempPb tempQb tempPc tempQc

%% Voltage angle approximation equality constraints

% If no phase present at node - No constraint on phase angle
% If phase present at node - Propagate voltage up feeder (toward head)
% |V_m||V_n|(-Theta_m + Theta_n) + 1/2 N_mn P_n - 1/2 M_mn Q_n = 0
% Above is not exact equation, read papers to find it

% Phase A
for k1 = 2:n
    if PH(1,k1) == 1
        tempDelta = -[FM(k1,:)==-1]; tempDelta(1,k1) = 1;
%         if sim.iter == 1
            j1 = find(FM(k1,:)==-1);
            tempDelta = tempDelta*abs(Vmagfbs(1,j1))*abs(Vmagfbs(1,k1));
%         end      
        tempPa = zeros(1,n); tempPa(1,k1) = -imag(gammaFZpu(1,1,k1));
        tempPb = zeros(1,n); tempPb(1,k1) = -imag(gammaFZpu(1,2,k1));
        tempPc = zeros(1,n); tempPc(1,k1) = -imag(gammaFZpu(1,3,k1));
        tempQa = zeros(1,n); tempQa(1,k1) = -real(gammaFZpu(1,1,k1));
        tempQb = zeros(1,n); tempQb(1,k1) = -real(gammaFZpu(1,2,k1));
        tempQc = zeros(1,n); tempQc(1,k1) = -real(gammaFZpu(1,3,k1));

        Aeq = [Aeq;
            zeros(1,n) zeros(1,n) pi/180*tempDelta tempPa tempQa zeros(1,n) zeros(1,n), ...
            zeros(1,n) zeros(1,n) zeros(1,n) tempPb tempQb zeros(1,n) zeros(1,n), ...
            zeros(1,n) zeros(1,n) zeros(1,n) tempPc tempQc zeros(1,n) zeros(1,n)];
        beq = [beq;
            0];
    end        
end

% Phase B
for k1 = 2:n
    if PH(2,k1) == 1
        tempDelta = -[FM(k1,:)==-1]; tempDelta(1,k1) = 1;
%         if sim.iter == 1
            j1 = find(FM(k1,:)==-1);
            tempDelta = tempDelta*abs(Vmagfbs(2,j1))*abs(Vmagfbs(2,k1));
%         end        
        tempPa = zeros(1,n); tempPa(1,k1) = -imag(gammaFZpu(2,1,k1));
        tempPb = zeros(1,n); tempPb(1,k1) = -imag(gammaFZpu(2,2,k1));
        tempPc = zeros(1,n); tempPc(1,k1) = -imag(gammaFZpu(2,3,k1));
        tempQa = zeros(1,n); tempQa(1,k1) = -real(gammaFZpu(2,1,k1));
        tempQb = zeros(1,n); tempQb(1,k1) = -real(gammaFZpu(2,2,k1));
        tempQc = zeros(1,n); tempQc(1,k1) = -real(gammaFZpu(2,3,k1));

        Aeq = [Aeq;
            zeros(1,n) zeros(1,n) zeros(1,n) tempPa tempQa zeros(1,n) zeros(1,n), ...
            zeros(1,n) zeros(1,n) pi/180*tempDelta tempPb tempQb zeros(1,n) zeros(1,n), ...
            zeros(1,n) zeros(1,n) zeros(1,n) tempPc tempQc zeros(1,n) zeros(1,n)];
        beq = [beq;
            0];
    end
end

% Phase C
for k1 = 2:n
    if PH(3,k1) == 1
        tempDelta = -[FM(k1,:)==-1]; tempDelta(1,k1) = 1;
%         if sim.iter == 1
            j1 = find(FM(k1,:)==-1);
            tempDelta = tempDelta*abs(Vmagfbs(3,j1))*abs(Vmagfbs(3,k1));
%         end        
        tempPa = zeros(1,n); tempPa(1,k1) = -imag(gammaFZpu(3,1,k1));
        tempPb = zeros(1,n); tempPb(1,k1) = -imag(gammaFZpu(3,2,k1));
        tempPc = zeros(1,n); tempPc(1,k1) = -imag(gammaFZpu(3,3,k1));
        tempQa = zeros(1,n); tempQa(1,k1) = -real(gammaFZpu(3,1,k1));
        tempQb = zeros(1,n); tempQb(1,k1) = -real(gammaFZpu(3,2,k1));
        tempQc = zeros(1,n); tempQc(1,k1) = -real(gammaFZpu(3,3,k1));

        Aeq = [Aeq;
            zeros(1,n) zeros(1,n) zeros(1,n) tempPa tempQa zeros(1,n) zeros(1,n), ...
            zeros(1,n) zeros(1,n) zeros(1,n) tempPb tempQb zeros(1,n) zeros(1,n), ...
            zeros(1,n) zeros(1,n) pi/180*tempDelta tempPc tempQc zeros(1,n) zeros(1,n)];
        beq = [beq;
            0];
    end
end

clear tempDelta tempPa tempQa tempPb tempQb tempPc tempQc

%% Power flow equality constraints

% If no phase present at node - Set power entering nodes to zero -
% P_n^phi = 0, Q_n^phi = 0
% If phase present at node -
% P_m^phi - u_m^phi - sum_n ( P_n^phi ) = p_n^phi + sum_n Re{L_mn^phi}
% Q_m^phi - v_m^phi - sum_n ( Q_n^phi ) = q_n^phi - cap_n^phi + sum_n Im{L_mn^phi}
% p_n^phi = p_n^phi*(A_PQ,n^phi + A_I,n^phi |V_n^phi| + A_Z,n^phi y_n^phi) + u_n^phi
% q_n^phi = q_n^phi*(A_PQ,n^phi + A_I,n^phi |V_n^phi| + A_Z,n^phi y_n^phi) -cap_n^phi + v_n^phi

for ph = 1:3
    for k1 = 1:n
        if PH(ph,k1) == 0
            tempPQ = zeros(1,n); tempPQ(1,k1) = 1;
            Aeq = [Aeq;
                    zeros(1,(ph-1)*nvar), zeros(1,3*n) tempPQ zeros(1,3*n), zeros(1,(3-ph)*nvar);
                    zeros(1,(ph-1)*nvar), zeros(1,4*n) tempPQ zeros(1,2*n), zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0;
                0];
        elseif PH(ph,k1) == 1
            tempY = zeros(1,n); tempY(1,k1) = -spu(ph,k1)*AZ(ph,k1);
            tempVmag = zeros(1,n); tempVmag(1,k1) = -spu(ph,k1)*AI(ph,k1);
            tempPQ = -[FM(k1,:)==1]; tempPQ(1,k1) = 1;
            tempuv = zeros(1,n); tempuv(1,k1) = -1;
            Aeq  = [Aeq;
                zeros(1,(ph-1)*nvar), real(tempY) real(tempVmag) zeros(1,n) tempPQ zeros(1,n) tempuv zeros(1,n), zeros(1,(3-ph)*nvar);
                zeros(1,(ph-1)*nvar), imag(tempY) imag(tempVmag) zeros(1,n) zeros(1,1*n) tempPQ zeros(1,n) tempuv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                real(spu(ph,k1)*APQ(ph,k1));
                imag(spu(ph,k1)*APQ(ph,k1)) - imag(cappu(ph,k1))];
            for k2 = 1:n
                if feeder.FM(k1,k2) == 1
                    beq(end-1:end) = beq(end-1:end) + [real(Lfbs(ph,k2)); imag(Lfbs(ph,k2))];
                end
            end
        end        
    end        
end

clear tempY tempPQ tempuv

%% Zero inverter outputs for non control nodes and nonexistent phases

% For all phases at all nodes that are non control, set inverter output to
% zero - u_n^phi = 0, v_n^phi = 0

for ph = 1:3
    for k1 = 1:n
        if cn(ph,k1) == 0 || PH(ph,k1) == 0 || cnstate(ph,k1) == 0
            tempuv = zeros(1,n); tempuv(1,k1) = 1;
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), zeros(1,5*n) tempuv zeros(1,n), zeros(1,(3-ph)*nvar);
                zeros(1,(ph-1)*nvar), zeros(1,5*n) zeros(1,n) tempuv, zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0;
                0];
        end
    end
end

clear tempuv

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

Acirc = zeros(2,3*nvar,3*n);
bcirc = zeros(3*n);

for ph = 1:3
    for k1 = 1:n
        Acirc(1, (ph-1)*nvar + 5*n + k1,(ph-1)*n + k1) = 1;
        Acirc(2, (ph-1)*nvar + 6*n + k1,(ph-1)*n + k1) = 1;
        if cn(ph,k1) == 1 && PH(ph,k1) == 0
            bcirc((ph-1)*n + k1) = 0;
        elseif  cn(ph,k1) == 1 && PH(ph,k1) == 1
            bcirc((ph-1)*n + k1) = wmaxpu(ph,k1);  
        end
    end
end

%% Inverter control bounds - hyperplane approximation of circle

% If no phase present at node - Set magnitude of inverter output to zero -
% |w_n^phi| <= 0
% If phase present at node - Set magnitude of inverter output to maximum
% for particular phase and node - |w_n^phi| <= wmax_n^phi

ncirc = 120;
theta = linspace(0,360,ncirc+1);

Aw = [];
bw = [];

for ph = 1:3
    for k1 = 2:n
        if controllers.cn(ph,k1) == 1
            tempu = zeros(ncirc+1,n); tempu(:,k1) = cosd(theta');
            tempv = zeros(ncirc+1,n); tempv(:,k1) = sind(theta');
            Aw = [Aw;
                zeros(ncirc+1,(ph-1)*nvar), zeros(ncirc+1,5*n) tempu tempv, zeros(ncirc+1,(3-ph)*nvar)];
            bw = [bw;
                wmaxpu(ph,k1)*ones(length(theta),1)];
        end
    end
end

%% CVX

tic

disp('CVX')

fy = zeros(3*nvar,1);
fy(nvar*(sim.tp-1)+sim.tn) = 1;
fy;

% Aineq = [Aineq;
%     Aw];
% bineq = [bineq;
%     bw];

cvx_begin quiet
    expressions Zy;
    variable X(3*nvar);
    Zy = fy'*X;
    minimize(Zy)
    subject to
    Aeq * X == beq
%     Aineq * X <= bineq
    [Aineq; Aw] * X <= [bineq; bw]
%     for ph = 1:3
%         for k1 = 1:n
%             if cn(ph,k1) == 1
%                 norm(Acirc(:,:,(ph-1)*n + k1)*X,2) <= bcirc((ph-1)*n + k1)
% %                 if sim.ramp == 1
% %                     norm(Aramp(:,:,(ph-1)*n + k1)*X - [real(wk1(ph,k1)); imag(wk1(ph,k1))],2) <= bramp((ph-1)*n + k1)
% %                 end
%             end
%         end
%     end
cvx_end

OPTSOLMIN.X = X;
OPTSOLMIN.Zy = Zy;
OPTSOLMIN.nvar = nvar;
OPTSOLMIN.cvx_status = cvx_status;

cvx_begin quiet
    expressions Zy;
    variable X(3*nvar);
    Zy = fy'*X;
    maximize(Zy)
    subject to
    Aeq * X == beq
%     Aineq * X <= bineq
    [Aineq; Aw] * X <= [bineq; bw]
%     for ph = 1:3
%         for k1 = 1:n
%             if cn(ph,k1) == 1
%                 norm(Acirc(:,:,(ph-1)*n + k1)*X,2) <= bcirc((ph-1)*n + k1)
% %                 if sim.ramp == 1
% %                     norm(Aramp(:,:,(ph-1)*n + k1)*X - [real(wk1(ph,k1)); imag(wk1(ph,k1))],2) <= bramp((ph-1)*n + k1)
% %                 end
%             end
%         end
%     end
cvx_end

OPTSOLMAX.X = X;
OPTSOLMAX.Zy = Zy;
OPTSOLMAX.nvar = nvar;
OPTSOLMAX.cvx_status = cvx_status;

toc

disp('LP MAKER')

tic


A = [Aineq; Aw; Aeq];
b = [bineq; bw; beq];
e = [-ones(size(Aineq,1),1); -ones(size(Aw,1),1); zeros(size(Aeq,1),1)];
vlb = -Inf*ones(length(fy),1);
vub = Inf*ones(length(fy),1);
vint = [];

lpmin = lp_maker(fy,A,b,e,vlb,vub,vint,1,1);
mxlpsol = mxlpsolve('solve',lpmin);
Zy = mxlpsolve('get_objective',lpmin);
X = mxlpsolve('get_variables',lpmin);
mxlpsolve('print_lp',lpmin);
mxlpsolve('delete_lp',lpmin);

OPTSOLMIN.X = X;
OPTSOLMIN.Zy = Zy;

lpmax = lp_maker(fy,A,b,e,vlb,vub,vint,1,0);
mxlpsol = mxlpsolve('solve',lpmax);
Zy = mxlpsolve('get_objective',lpmax);
X = mxlpsolve('get_variables',lpmax);
mxlpsolve('print_lp',lpmax);
mxlpsolve('delete_lp',lpmax);

OPTSOLMAX.X = X;
OPTSOLMAX.Zy = Zy;


toc

end
