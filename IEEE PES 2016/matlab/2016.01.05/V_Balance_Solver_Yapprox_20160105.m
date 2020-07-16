function [X,Y,nvar,cvx_optval,cvx_status] = V_Balance_Solver_Yapprox_20160105(feeder,loads,controllers)
%%

% Michael Sankur - msankur@berkeley.edu
% 2016.01.05

% Variables - [Y_a; P_a; Q_a; u_a; v_a;
%               Y_b; P_b; Q_b; u_b; v_b;
%               Y_c; P_c; Q_c; u_c; v_c];

%%

% feeder parameters
n = feeder.n;
FM = feeder.FM;
PH = feeder.PH;
FZpu = feeder.FZpu;

% load parameters
spu = loads.spu;
A0 = loads.A0;
A1 = loads.A1;
cappu = loads.cappu;

% controller parameters
cn = controllers.cn;
cnstate = controllers.cnstate;
wmaxpu = controllers.wmaxpu;
rrpu = controllers.rrpu;

%%

nvar = 5*n; % number of variables per phase

%% Set up optimization matrices

Aineq = []; % Inequality constraint matrix
bineq = []; % Inequality constraint vector
Aeq = []; % Equality constraint matrix
beq = []; % Equality constraint vector

%% Voltage magnitude equality and inequality constraints

% If no phase present at node - Set voltage magnitude to zero
% If phase present at node - Add inequality constraint on squared voltage
% magnitude - 0.95^2 <= y_m,k <= 1.05^2

for ph = 1:3
    for k1 = 2:n
        if PH(ph,k1) == 0
            tempY = zeros(1,n); tempY(1,k1) = 1;
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), tempY zeros(1,4*n), zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0];
        elseif PH(ph,k1) == 1
            tempY = zeros(2,n); tempY(:,k1) = [1; -1];
            Aineq = [Aineq;
                zeros(2,(ph-1)*nvar), tempY zeros(2,4*n), zeros(2,(3-ph)*nvar)];
            bineq = [bineq;
                1.05^2;
                -(0.95^2)];
        end
    end
end

clear tempY

%% Feeder head voltage magnitude (initialize) constraints

% Initialize feeder head square voltage magnitude as 1 for all three phases
% y_m,1 = 1

tempY = zeros(1,n); tempY(1) = 1;

Aeq = [Aeq;
    tempY zeros(1,4*n), zeros(1,nvar) zeros(1,nvar);
    zeros(1,nvar), tempY zeros(1,4*n), zeros(1,nvar);
    zeros(1,nvar), zeros(1,nvar), tempY zeros(1,4*n)];

beq = [beq;
    1;
    1;
    1];

clear tempY

%% Voltage magnitude approximation equality constraints

% If no phase present at node - Set voltage at phase to zero - y_m,k = 0
% If phase present at node - Propagate voltage up feeder (toward head)
% -Y_j + Y_k + M_jk^P P_k + M_jk^Q Q_k = 0

% Phase A
for k1 = 2:n
    if PH(1,k1) == 0
%         tempY = zeros(1,n); tempY(1,k1) = 1;
%         Aeq = [Aeq;
%             tempY zeros(1,4*n), zeros(1,nvar), zeros(1,nvar)];
%         beq = [beq;
%             0];       
    elseif PH(1,k1) == 1
        tempY = -[FM(k1,:)==-1]; tempY(1,k1) = 1;
        tempPa = zeros(1,n); tempPa(1,k1) = 2*real(FZpu(1,1,k1));
        tempPb = zeros(1,n); tempPb(1,k1) = -real(FZpu(1,2,k1)) + sqrt(3)*imag(FZpu(1,2,k1));
        tempPc = zeros(1,n); tempPc(1,k1) = -real(FZpu(1,3,k1)) - sqrt(3)*imag(FZpu(1,3,k1));
        tempQa = zeros(1,n); tempQa(1,k1) = 2*imag(FZpu(1,1,k1));    
        tempQb = zeros(1,n); tempQb(1,k1) = -sqrt(3)*real(FZpu(1,2,k1)) - imag(FZpu(1,2,k1)) ;    
        tempQc = zeros(1,n); tempQc(1,k1) = sqrt(3)*real(FZpu(1,3,k1)) - imag(FZpu(1,3,k1));

        Aeq = [Aeq;
            tempY tempPa tempQa zeros(1,2*n), ...
            zeros(1,n) tempPb tempQb zeros(1,2*n), ...
            zeros(1,n) tempPc tempQc zeros(1,2*n)];
        beq = [beq;
            0];
    end      
end

% Phase B
for k1 = 2:n
    if PH(2,k1) == 0
%         tempY = zeros(1,n); tempY(1,k1) = 1;
%         Aeq = [Aeq;
%             zeros(1,nvar), tempY zeros(1,4*n), zeros(1,nvar)];
%         beq = [beq;
%             0];
    elseif PH(2,k1) == 1
        tempY = -[FM(k1,:)==-1]; tempY(1,k1) = 1;
        tempPa = zeros(1,n); tempPa(1,k1) = -real(FZpu(2,1,k1)) - sqrt(3)*imag(FZpu(2,1,k1));
        tempPb = zeros(1,n); tempPb(1,k1) = 2*real(FZpu(2,2,k1));
        tempPc = zeros(1,n); tempPc(1,k1) = -real(FZpu(2,3,k1)) + sqrt(3)*imag(FZpu(2,3,k1));
        tempQa = zeros(1,n); tempQa(1,k1) = sqrt(3)*real(FZpu(2,1,k1)) - imag(FZpu(2,1,k1));    
        tempQb = zeros(1,n); tempQb(1,k1) = 2*imag(FZpu(2,2,k1));    
        tempQc = zeros(1,n); tempQc(1,k1) = -sqrt(3)*real(FZpu(2,3,k1)) - imag(FZpu(2,3,k1));

        Aeq = [Aeq;
            zeros(1,n) tempPa tempQa zeros(1,2*n), ...
            tempY tempPb tempQb zeros(1,2*n), ...
            zeros(1,n) tempPc tempQc zeros(1,2*n)];
        beq = [beq;
            0];
    end
end

% Phase C
for k1 = 2:n
    if PH(3,k1) ==0
%         tempY = zeros(1,n); tempY(1,k1) = 1;
%         Aeq = [Aeq;
%             zeros(1,nvar), zeros(1,nvar), tempY zeros(1,4*n)];
%         beq = [beq;
%             0];
    elseif PH(3,k1) == 1
        tempY = -[FM(k1,:)==-1]; tempY(1,k1) = 1;
        tempPa = zeros(1,n); tempPa(1,k1) = -real(FZpu(3,1,k1)) + sqrt(3)*imag(FZpu(3,1,k1));
        tempPb = zeros(1,n); tempPb(1,k1) = -real(FZpu(3,2,k1)) - sqrt(3)*imag(FZpu(3,2,k1));
        tempPc = zeros(1,n); tempPc(1,k1) = 2*real(FZpu(3,3,k1));
        tempQa = zeros(1,n); tempQa(1,k1) = -sqrt(3)*real(FZpu(3,1,k1)) - imag(FZpu(3,1,k1));    
        tempQb = zeros(1,n); tempQb(1,k1) = sqrt(3)*real(FZpu(3,2,k1)) - imag(FZpu(3,2,k1));
        tempQc = zeros(1,n); tempQc(1,k1) = 2*imag(FZpu(3,3,k1));

        Aeq = [Aeq;
            zeros(1,n) tempPa tempQa zeros(1,2*n), ...
            zeros(1,n) tempPb tempQb zeros(1,2*n), ...
            tempY tempPc tempQc zeros(1,2*n)];
        beq = [beq;
            0];
    end
end

clear tempY tempPa tempQa tempPb tempQb tempPc tempQc

%% Power flow equality constraints

% If no phase present at node - Set power entering nodes to zero -
% P_m,k = 0, Q_m,k = 0
% If phase present at node -
% If non-control node - Conservation of power, ignoring inverter output
% P_m,k - u_m,k - sum_j ( P_j,k ) = p_m,k
% Q_m,k - v_m,k - sum_j ( Q_j,k ) = q_m,k - cap_m,k

for ph = 1:3
    for k1 = 1:n
        if PH(ph,k1) == 0
            tempPQ = zeros(1,n); tempPQ(1,k1) = 1;
            Aeq = [Aeq;
                    zeros(1,(ph-1)*nvar), zeros(1,1*n) tempPQ zeros(1,3*n), zeros(1,(3-ph)*nvar);
                    zeros(1,(ph-1)*nvar), zeros(1,2*n) tempPQ zeros(1,2*n), zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0;
                0];
        elseif PH(ph,k1) == 1
%             % Non control nodes
%             if cn(k1) == 0
%                 tempY = zeros(1,n); tempY(1,k1) = spu(ph,k1)*A1(ph,k1);
%                 tempPQ = -[FM(k1,:)==1]; tempPQ(1,k1) = 1;
%                 Aeq = [Aeq;
%                     zeros(1,(ph-1)*nvar), real(tempY) tempPQ zeros(1,3*n), zeros(1,(3-ph)*nvar);
%                     zeros(1,(ph-1)*nvar), imag(tempY) zeros(1,1*n) tempPQ zeros(1,2*n), zeros(1,(3-ph)*nvar)];
%             % Control nodes
%             elseif cn(k1) == 1
                tempY = zeros(1,n); tempY(1,k1) = -spu(ph,k1)*A1(ph,k1);
                tempPQ = -[FM(k1,:)==1]; tempPQ(1,k1) = 1;
                tempuv = zeros(1,n); tempuv(1,k1) = -1;
                Aeq  = [Aeq;
                    zeros(1,(ph-1)*nvar), real(tempY) tempPQ zeros(1,n) tempuv zeros(1,n), zeros(1,(3-ph)*nvar);
                    zeros(1,(ph-1)*nvar), imag(tempY) zeros(1,1*n) tempPQ zeros(1,n) tempuv, zeros(1,(3-ph)*nvar)];
%             end
            beq = [beq;
                real(spu(ph,k1)*A0(ph,k1));
                imag(spu(ph,k1)*A0(ph,k1)) - imag(cappu(ph,k1))];
        end        
    end        
end

clear tempY tempPQ tempuv

%% Zero inverter outputs for non control nodes and nonexistent phases

% For all phases at all nodes that are non control, set inverter output to
% zero - u_m,k = 0, v_m,k = 0

for ph = 1:3
    for k1 = 1:n
        if cn(ph,k1) == 0 || PH(ph,k1) == 0 || cnstate(ph,k1) == 0
            tempuv = zeros(1,n); tempuv(1,k1) = 1;
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), zeros(1,3*n) tempuv zeros(1,n), zeros(1,(3-ph)*nvar);
                zeros(1,(ph-1)*nvar), zeros(1,3*n) zeros(1,n) tempuv, zeros(1,(3-ph)*nvar)];
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

%% Inverter ability to consume real power constraints

% for ph = 1:3
%     for k1 = 1:n
%         if cn(k1) == 1
%             if PH(ph,k1) == 0
% 
%             elseif PH(ph,k1) == 1
%                 tempuv = zeros(1,n); tempuv(:,k1) = 1;
%                 Aineq = [Aineq;
%                     zeros(1,(ph-1)*nvar), zeros(1,3*n) tempuv zeros(1,n), zeros(1,(3-ph)*nvar);
%                 bineq = [bineq;
%                     0];
%             end                       
%         end
%     end
% end
% 
% clear tempuv

%% Inverter ability to only supply reactive power constraints

% u_m,k = 0

for ph = 1:3
    for k1 = 1:n
        if cn(ph,k1) == 1 && PH(ph,k1) == 1 && cnstate(ph,k1) == 1
            tempuv = zeros(1,n); tempuv(:,k1) = 1;
            Aeq = [Aeq;
                zeros(1,(ph-1)*nvar), zeros(1,3*n) tempuv zeros(1,n), zeros(1,(3-ph)*nvar)];
            beq = [beq;
                0];                    
        end
    end
end

clear tempuv

%% Inverter reactive power constraints

% -wmax(k) <= v_m,k <= wmaxpu(k)

% for ph = 1:3
%     for k1 = 1:n
%         if cn(k1) == 1 && PH(ph,k1) == 1
%             tempuv = zeros(2,n); tempuv(:,k1) = [1; -1];
%             Aineq = [Aineq;
%                 zeros(2,(ph-1)*nvar), zeros(2,4*n) tempuv zeros(2,0*n), zeros(2,(3-ph)*nvar)];
%             bineq = [bineq;
%                 wmaxpu(k1);
% %                     0;
%                 wmaxpu(k1)];                  
%         end
%     end
% end
% 
% clear tempuv


%% Inverter control bounds - circle

% If no phase present at node - Set magnitude of inverter output to zero -
% |w_m,k| <= 0
% If phase present at node - Set magnitude of inverter output to maximum
% for particular phase and node - |w_m,k| <= w_max(m,k)

Acirc = zeros(3*nvar,3*nvar,3*n);
bcirc = zeros(n);

for ph = 1:3
    for k1 = 1:n
        Acirc((ph-1)*nvar + 3*n + k1, (ph-1)*nvar + 3*n + k1,(ph-1)*n + k1) = 1;
        Acirc((ph-1)*nvar + 4*n + k1, (ph-1)*nvar + 4*n + k1,(ph-1)*n + k1) = 1;
        if cn(ph,k1) == 1 && PH(ph,k1) == 0
            bcirc((ph-1)*n + k1) = 0;
        elseif  cn(ph,k1) == 1 && PH(ph,k1) == 1
            bcirc((ph-1)*n + k1) = wmaxpu(ph,k1);  
        end
    end
end

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
        FV(:,:,k1) = [tempY zeros(1,4*n), -tempY zeros(1,4*n), zeros(1,5*n);
            zeros(2,nvar), zeros(2,nvar), zeros(2,nvar)];
    elseif PH(1,k1) == 0 && PH(2,k1) == 1 && PH(3,k1) == 1
        FV(:,:,k1) = [zeros(1,5*n), tempY zeros(1,4*n), -tempY zeros(1,4*n);
            zeros(2,nvar), zeros(2,nvar), zeros(2,nvar)];
    elseif PH(1,k1) == 1 && PH(2,k1) == 0 && PH(3,k1) == 1
        FV(:,:,k1) = [tempY zeros(1,4*n), zeros(1,5*n), -tempY zeros(1,4*n);
            zeros(2,nvar), zeros(2,nvar), zeros(2,nvar)];   
    elseif PH(1,k1) == 1 && PH(2,k1) == 1 && PH(3,k1) == 1
        FV(:,:,k1) = [tempY zeros(1,4*n), -tempY zeros(1,4*n), zeros(1,nvar);
            zeros(1,nvar), tempY zeros(1,4*n), -tempY zeros(1,4*n);
            tempY zeros(1,4*n), zeros(1,nvar), -tempY zeros(1,4*n)];
    end
end
FV;
size(FV);

% Matrix weighting control effort
Fu1 = blkdiag(zeros(3*n,3*n),diag(controllers.cn(1,:)),diag(controllers.cn(1,:)));
Fu2 = blkdiag(zeros(3*n,3*n),diag(controllers.cn(2,:)),diag(controllers.cn(2,:)));
Fu3 = blkdiag(zeros(3*n,3*n),diag(controllers.cn(3,:)),diag(controllers.cn(3,:)));
Fu = blkdiag(Fu1,Fu2,Fu3);

% Vector for feeder head power in objective function
fP0 = zeros(3*nvar,1);
for ph = 1:3
    fP0(1*n + 1 + (ph-1)*nvar,1) = 1;
end
fP0;


cvx_begin quiet;
    expression Z;
    variable X(3*nvar);
    dual variable Y;
    for k1 = 2:13
        Z = Z +  norm(FV(:,:,k1)*X,2);
    end
    Z = Z + 0.5*norm(Fu*X,2);
    minimize(Z)
    subject to;
    Aeq * X == beq;
    Y : Aineq * X <= bineq;
    for ph = 1:3
        for k1 = 1:n
            if cn(ph,k1) == 1
                norm(Acirc(:,:,(ph-1)*n + k1)*X,2) <= bcirc((ph-1)*n + k1);
            end
        end
    end
cvx_end;


end
