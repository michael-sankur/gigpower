function OPTSOL = Solver_LinDist3Flow_Tracking_20160603(feeder,loads,controllers,sim)
%%

% Michael Sankur - msankur@berkeley.edu, msankur@lbl.gov
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

mag_ref = sim.mag_ref;
delta_ref = sim.delta_ref;

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
    gammaFZpu(:,:,k1) = gamma(:,:,k1).*(FZpu(:,:,k1)');
end
gamma;
gammaFZpu;

%% Voltage magnitude approximation equality constraints

% If no phase present at node - Set voltage at phase to zero - y_n^phi = 0
% If phase present at node - Propagate voltage up feeder (toward head)
% -Y_m + Y_n + M_mn P_n + N_mn Q_n = -H_mn

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
                zeros(1,(ph-1)*nvar), imag(tempY) real(tempVmag) zeros(1,n) zeros(1,1*n) tempPQ zeros(1,n) tempuv, zeros(1,(3-ph)*nvar)];
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
bcirc = zeros(n);

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

Aramp = zeros(2,3*nvar,3*n);
bramp = zeros(3*n,1);

for ph = 1:3
    for k1 = 1:n
        Aramp(1,(ph-1)*nvar + 5*n + k1,(ph-1)*n + k1) = 1;
        Aramp(2,(ph-1)*nvar + 6*n + k1,(ph-1)*n + k1) = 1;
        if cn(ph,k1) == 1 && PH(ph,k1) == 0
            bramp((ph-1)*n + k1) = 0;
        elseif cn(ph,k1) == 1 && PH(ph,k1) == 1
            bramp((ph-1)*n + k1) = rrpu(ph,k1);
        end
    end
end

clear tempuv

%% Constraint for magnitude

% for ph = 1:3
%     for k1 = 1:n
%         if isnan(delta_ref(ph,k1)) == 0
%             templine = zeros(1,3*nvar);
%             templine(1,(ph-1)*nvar + k1) = 1;
%             Aeq = [Aeq;
%                 templine];
%             beq = [beq;
%                 mag_ref(ph,k1)^2];
%         end
%     end
% end

%% Constraint for phase angle

% for ph = 1:3
%     for k1 = 1:n
%         if isnan(delta_ref(ph,k1)) == 0
%             templine = zeros(1,3*nvar);
%             templine(1,(ph-1)*nvar + 2*n + k1) = 1;
%             Aeq = [Aeq;
%                 templine];
%             beq = [beq;
%                 delta_ref(ph,k1)];
%         end
%     end
% end

%% Objective function - corrected

% Matrices to create a scalar or vector of voltage magnitude differences
% for each node. If a node has one phase, the vector is null. If the node
% has two phases, the vector contains Y_m,k - Y_n,k where m and n are the
% existing phases. If the node has three phases, the vector is as follows:
% [Y_a,k - Y_b,k ; Y_b,k - Y_c,k ; Y_a,k - Y_c,k].
FV = [];
for k1 = 2:n
    tempY = zeros(1,n); tempY(k1) = 1;
    if sum(PH(:,k1)) == 1
        FV(:,:,k1) = zeros(3,3*nvar);
    elseif PH(1,k1) == 1 && PH(2,k1) == 1 && PH(3,k1) == 0
        FV(:,:,k1) = [tempY zeros(1,6*n), -tempY zeros(1,6*n), zeros(1,nvar);
            zeros(2,nvar), zeros(2,nvar), zeros(2,nvar)];
    elseif PH(1,k1) == 0 && PH(2,k1) == 1 && PH(3,k1) == 1
        FV(:,:,k1) = [zeros(1,nvar), tempY zeros(1,6*n), -tempY zeros(1,6*n);
            zeros(2,nvar), zeros(2,nvar), zeros(2,nvar)];
    elseif PH(1,k1) == 1 && PH(2,k1) == 0 && PH(3,k1) == 1
        FV(:,:,k1) = [tempY zeros(1,6*n), zeros(1,nvar), -tempY zeros(1,6*n);
            zeros(2,nvar), zeros(2,nvar), zeros(2,nvar)];   
    elseif PH(1,k1) == 1 && PH(2,k1) == 1 && PH(3,k1) == 1
        FV(:,:,k1) = [tempY zeros(1,6*n), -tempY zeros(1,6*n), zeros(1,nvar);
            zeros(1,nvar), tempY zeros(1,6*n), -tempY zeros(1,6*n);
            tempY zeros(1,6*n), zeros(1,nvar), -tempY zeros(1,6*n)];
    end
end
FV;
size(FV);

% Matrix weighting control effort
Fu1 = blkdiag(zeros(5*n,5*n),diag(controllers.cn(1,:)),diag(controllers.cn(1,:)));
Fu2 = blkdiag(zeros(5*n,5*n),diag(controllers.cn(2,:)),diag(controllers.cn(2,:)));
Fu3 = blkdiag(zeros(5*n,5*n),diag(controllers.cn(3,:)),diag(controllers.cn(3,:)));
Fu = blkdiag(Fu1,Fu2,Fu3);

% Matrix weigting control effort
FU = [];
for k1 = 2:n
    tempuv = zeros(1,n); tempuv(k1) = 1;
    tempphase = [zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) tempuv zeros(1,n);
        zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) zeros(1,n) tempuv];
    FU(:,:,k1) = [tempphase zeros(2,nvar) zeros(2,nvar);
        zeros(2,nvar) tempphase zeros(2,nvar);
        zeros(2,nvar) zeros(2,nvar) tempphase];
end

% Vector for feeder head power in objective function
fP0 = zeros(nvar,1);
fP0(3*n+1,1) = 1;
fP0 = [fP0; fP0; fP0];

%% CVX

rho_y = 1e3;
rho_theta = 1e2;
rho_w = 1e0;

cvx_begin quiet
    expressions Zy Ztheta Zw;
    variable X(3*nvar);
    for ph = 1:3
        for k1 = 2:n
            if isnan(mag_ref(ph,k1)) == 0
                Zy = Zy + (X((ph-1)*nvar+k1,1) - mag_ref(ph,k1)^2)^2;
            end
            if isnan(delta_ref(ph,k1)) == 0
                Ztheta = Ztheta + (X((ph-1)*nvar+2*n+k1,1) - delta_ref(ph,k1))^2;
            end
            if controllers.cn(ph,k1) == 1
%                 Zw = Zw + rho_w*norm(FU(ph*2-1:ph*2,:,k1)*X,2);
                Zw = Zw + (X((ph-1)*nvar+5*n+k1)^2 + X((ph-1)*nvar+6*n+k1)^2);
            end
        end
    end
    minimize(rho_y*Zy + rho_theta*Ztheta + rho_w*Zw)
    subject to
    Aeq * X == beq
    Aineq * X <= bineq
    for ph = 1:3
        for k1 = 1:n
            if cn(ph,k1) == 1
                norm(Acirc(:,:,(ph-1)*n + k1)*X,2) <= bcirc((ph-1)*n + k1)
%                 if sim.ramp == 1
%                     norm(Aramp(:,:,(ph-1)*n + k1)*X - [real(wk1(ph,k1)); imag(wk1(ph,k1))],2) <= bramp((ph-1)*n + k1)
%                 end
            end
        end
    end
cvx_end

OPTSOL.X = X;
% OPTSOL.Y = Y;
OPTSOL.nvar = nvar;
OPTSOL.cvx_optval = cvx_optval;
OPTSOL.cvx_status = cvx_status;
OPTSOL.Zy = rho_y*Zy;
OPTSOL.Ztheta = rho_theta*Ztheta;
OPTSOL.Zw = rho_w*Zw;

end
