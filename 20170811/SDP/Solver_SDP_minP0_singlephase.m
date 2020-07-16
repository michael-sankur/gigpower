function [sdpsol,SDPMAT] = Solver_SDP_minP0_singlephase(feeder, nodes, lines, configs, loads, caps, controllers, sim)
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
NPH = nodes.PH;

% line paramters
nline = lines.nline;
TXnum = lines.TXnum;
RXnum = lines.RXnum;
FZpu = lines.FZpu;
FYpu = lines.FYpu;

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

slacknode = sim.slacknode;
V0 = sim.Vnom;

rho = sim.rho

% iteration parameters
% if isempty(sim.Vfbs)
%     Vfbs = [1*ones(1,nnode);
%         1*exp(j*-120*pi/180)*ones(1,nnode);
%         1*exp(j*120*pi/180)*ones(1,nnode)].*nodes.PH;
% else
%     Vfbs = sim.Vfbs;
% end
% if isempty(sim.Lfbs)
%     Lfbs = zeros(3,nline);
% else
%     Lfbs = sim.Lfbs;
% end
% if isempty(sim.Hfbs)
%     Hfbs = zeros(3,nline);
% else
%     Hfbs1 = sim.Hfbs;
% end

nvar = nnode + nnode + nline + nline + nnode + nnode; % number of variables per phase

%% Convert to single phase


%% Nodes parameters

NPH = nodes.PH(1,:);

%% Line parameters

for k1 = 1:nline
    tempFZpu(k1) = FZpu(1,1,k1);
    tempFYpu(k1) = FYpu(1,1,k1);
end

FZpu = tempFZpu;
FYpu = tempFYpu;

%% Load parameters

aPQ = 1.0*ones(1,nnode).*NPH;
aI = 0.0*ones(1,nnode).*NPH;
aZ = 0.0*ones(1,nnode).*NPH;

loads.spu = 1.5*loads.spu(1,:);

%% Capacitor parameters

cappu = cappu(1,:);

%% Controller parameters

controllers.wmaxpu = 1.0*controllers.wmaxpu(1,:);

%% Simulation parameters

V0 = V0(1);


%% Construct Y matrix

Y = zeros(nnode,nnode);

for k1 = 1:nline
    
    txn = lines.TXnum(k1);
    rxn = lines.RXnum(k1);
    
    Y(txn,txn) = Y(txn,txn) - FYpu(k1);
    Y(rxn,rxn) = Y(rxn,rxn) - FYpu(k1);
    Y(txn,rxn) = FYpu(k1);
    Y(rxn,txn) = FYpu(k1);    
end

Y == Y.';

% max(max(abs(Y - Y)))

%% Construct e vectors

ea = eye(nnode,nnode);
% for k1 = 1:nnode
%     ea(:,k1) = e(:,3*k1-2);
%     eb(:,k1) = e(:,3*k1-1);
%     ec(:,k1) = e(:,3*k1);
% end
% ea, eb, ec

%% Construct Phi and Y matrices

for k1 = 1:nnode
    
    PhiVa(:,:,k1) = ea(:,k1)*ea(:,k1).';
%     PhiVb(:,:,k1) = eb(:,k1)*eb(:,k1).';
%     PhiVc(:,:,k1) = ec(:,k1)*ec(:,k1).';
    
    Ya(:,:,k1) = Y'*PhiVa(:,:,k1);
%     Yb(:,:,k1) = Y'*PhiVb(:,:,k1);
%     Yc(:,:,k1) = Y'*PhiVc(:,:,k1);
    
%     Ya(:,:,k1) = PhiVa(:,:,k1)*Y;
%     Yb(:,:,k1) = PhiVb(:,:,k1)*Y;
%     Yc(:,:,k1) = PhiVc(:,:,k1)*Y;
    
    PhiIa(:,:,k1) = Y'*ea(:,k1)*ea(:,k1).'*Y;
%     PhiIb(:,:,k1) = Y'*eb(:,k1)*eb(:,k1).'*Y;
%     PhiIc(:,:,k1) = Y'*ec(:,k1)*ec(:,k1).'*Y;
       
    PhiPa(:,:,k1) = 0.5*(Ya(:,:,k1) + Ya(:,:,k1)');
%     PhiPb(:,:,k1) = 0.5*(Yb(:,:,k1) + Yb(:,:,k1)');
%     PhiPc(:,:,k1) = 0.5*(Yc(:,:,k1) + Yc(:,:,k1)');
    
    PhiPa(:,:,k1) = PhiPa(:,:,k1) - real(spu(1,k1))*aZ(1,k1)*PhiVa(:,:,k1);
%     PhiPb(:,:,k1) = PhiPb(:,:,k1) - aZ(2,k1)*real(spu(2,k1))*PhiVb(:,:,k1);
%     PhiPc(:,:,k1) = PhiPc(:,:,k1) - aZ(3,k1)*real(spu(3,k1))*PhiVc(:,:,k1);
    
    PhiQa(:,:,k1) = -1j*0.5*(Ya(:,:,k1) - Ya(:,:,k1)');
%     PhiQb(:,:,k1) = -1j*0.5*(Yb(:,:,k1) - Yb(:,:,k1)');
%     PhiQc(:,:,k1) = -1j*0.5*(Yc(:,:,k1) - Yc(:,:,k1)');
    
    PhiQa(:,:,k1) = PhiQa(:,:,k1) - imag(spu(1,k1))*aZ(1,k1)*PhiVa(:,:,k1);
%     PhiQb(:,:,k1) = PhiQb(:,:,k1) - aZ(2,k1)*imag(spu(2,k1))*PhiVb(:,:,k1);
%     PhiQc(:,:,k1) = PhiQc(:,:,k1) - aZ(3,k1)*imag(spu(3,k1))*PhiVc(:,:,k1);
    
end

%% Construct Phase Angle Matrices

for k1 = 1:nnode
    
    Phim1(:,:,k1) = ea(:,1)*ea(:,k1).';
    
    PhiRe(:,:,k1) = 0.5*(Phim1(:,:,k1) + Phim1(:,:,k1).');
    
    PhiIm(:,:,k1) = (0.5/1j)*(Phim1(:,:,k1) - Phim1(:,:,k1).');
    
end

%%

SDPMAT.Y = Y;

SDPMAT.ea = ea;
% SDPMAT.eb = eb;
% SDPMAT.ec = ec;

SDPMAT.PhiVa = PhiVa;
% SDPMAT.PhiVb = PhiVb;
% SDPMAT.PhiVc = PhiVc;

SDPMAT.Ya = Ya;
% SDPMAT.Yb = Yb;
% SDPMAT.Yc = Yc;

SDPMAT.PhiPa = PhiPa;
% SDPMAT.PhiPb = PhiPb;
% SDPMAT.PhiPc = PhiPc;

SDPMAT.PhiQa = PhiQa;
% SDPMAT.PhiQb = PhiQb;
% SDPMAT.PhiQc = PhiQc;

% SDPMAT.PhiAnga = PhiAnga;
% SDPMAT.PhiAngb = PhiAngb;
% SDPMAT.PhiAngc = PhiAngc;

%%

wmaxpu(1,1) = 2*ones(1,1);

eps = 1.0

cvx_begin sdp
    cvx_solver sedumi;
    variable Xsdp(nnode,nnode) hermitian;
    variable gammasdp;
    variable alphasdp;
    variable gammaw(1,nnode);
%     variable gammasdp(1,1);
    expression Zsdp;
    expression Zsdp2;
    Zsdp = -trace(PhiVa(:,:,3)*Xsdp);
%     Zsdp = trace(PhiPa(:,:,1)*Xsdp) + 5*trace(PhiPa(:,:,4)*Xsdp) + 1*trace(PhiPa(:,:,5)*Xsdp);
    for k1 = 2:nnode
        Zsdp2 = Zsdp2 - eps*(trace(PhiQa(:,:,k1)*Xsdp) - aPQ(1,k1)*imag(spu(1,k1)) + cappu(1,k1));
    end
%     Zsdp = gammasdp;
    minimize(Zsdp + 0*Zsdp2);
%     minimize(alphasdp)
    subject to
       
    for k1 = 1:nnode

        if wmaxpu(1,k1) == 0

            trace(PhiPa(:,:,k1)*Xsdp) == aPQ(1,k1)*real(spu(1,k1));
            trace(PhiQa(:,:,k1)*Xsdp) == aPQ(1,k1)*imag(spu(1,k1)) - cappu(1,k1);

        else

%             trace(-PhiPa(:,:,k1)*Xsdp) >= -aPQ(1,k1)*real(spu(1,k1)) - wmaxpu(1,k1);
%             trace(PhiPa(:,:,k1)*Xsdp) >= aPQ(1,k1)*real(spu(1,k1)) - wmaxpu(1,k1);
% 
%             trace(-PhiQa(:,:,k1)*Xsdp) >= -aPQ(1,k1)*imag(spu(1,k1)) + cappu(1,k1) - wmaxpu(1,k1);
%             trace(PhiQa(:,:,k1)*Xsdp) >= aPQ(1,k1)*imag(spu(1,k1)) - cappu(1,k1) - wmaxpu(1,k1);

%         [-gammaw(1,k1) (trace(PhiPa(:,:,k1)*Xsdp) - aPQ(1,k1)*real(spu(1,k1))) (trace(PhiQa(:,:,k1)*Xsdp) - aPQ(1,k1)*imag(spu(1,k1)) + cappu(1,k1))
%             (trace(PhiPa(:,:,k1)*Xsdp) - aPQ(1,k1)*real(spu(1,k1))) -1 0;
%             (trace(PhiQa(:,:,k1)*Xsdp) - aPQ(1,k1)*imag(spu(1,k1)) + cappu(1,k1)) 0 -1] <= 0;

        [-wmaxpu(1,k1)^2 (trace(PhiPa(:,:,k1)*Xsdp) - aPQ(1,k1)*real(spu(1,k1))) (trace(PhiQa(:,:,k1)*Xsdp) - aPQ(1,k1)*imag(spu(1,k1)) + cappu(1,k1))
            (trace(PhiPa(:,:,k1)*Xsdp) - aPQ(1,k1)*real(spu(1,k1))) -1 0;
            (trace(PhiQa(:,:,k1)*Xsdp) - aPQ(1,k1)*imag(spu(1,k1)) + cappu(1,k1)) 0 -1] <= 0;

        end

        trace(-PhiVa(:,:,k1)*Xsdp) >= -1.05^2;
        trace(PhiVa(:,:,k1)*Xsdp) >= 0.95^2;
        
%         for delta = 0:5:360        
%             (trace(PhiPa(:,:,k1)*Xsdp) - real(spu(1,k1)))*cosd(delta) + ...
%                 (trace(PhiQa(:,:,k1)*Xsdp) - imag(spu(1,k1)) + cappu(1,k1))*sind(delta) <= wmaxpu(1,k1);
%         end
        
    end
    
%     k1 = 3;
%     Aw = [-wmaxpu(1,k1)^2 (trace(PhiPa(:,:,k1)*Xsdp) - aPQ(1,k1)*real(spu(1,k1))) (trace(PhiQa(:,:,k1)*Xsdp) - aPQ(1,k1)*imag(spu(1,k1)) + cappu(1,k1))
%         (trace(PhiPa(:,:,k1)*Xsdp) - aPQ(1,k1)*real(spu(1,k1))) -1 0;
%         (trace(PhiQa(:,:,k1)*Xsdp) - aPQ(1,k1)*imag(spu(1,k1)) + cappu(1,k1)) 0 -1];
%     Aw <= 0;
        
    
%     Magnitude Reference
    
%     trace(PhiVa(:,:,3)*Xsdp) - 0.965^2 <= gammasdp;
%     trace(-PhiVa(:,:,3)*Xsdp) + 0.965^2 <= gammasdp;

    % Magnitude Difference Minimization
    
%     trace(PhiVa(:,:,4)*Xsdp) - trace(PhiVa(:,:,6)*Xsdp) <= gammasdp;
%     -trace(PhiVa(:,:,4)*Xsdp) + trace(PhiVa(:,:,6)*Xsdp) <= gammasdp;
    
    % Angle Absolute Bounds
    
%     lb = -0.85;    
%     trace(PhiRe(:,:,4)*Xsdp)*tand(lb) <= trace(PhiIm(:,:,4)*Xsdp)
%     
%     ub = -0.65;    
%     trace(PhiIm(:,:,4)*Xsdp) <= trace(PhiRe(:,:,4)*Xsdp)*tand(ub)

    % Angle Difference Bounds
    
%     Phimn = ea(:,4)*ea(:,6).';
%     PhiRemn = 0.5*(Phimn + Phimn');
%     PhiImmn = (0.5/1j)*(Phimn - Phimn');
%     
%     lb = -0.1;
% %     trace(PhiRemn*Xsdp)*tand(lb) <= trace(PhiImmn*Xsdp)
%     Philb = PhiRemn*tand(lb) - PhiImmn;
%     trace(Philb*Xsdp) <= 0;
%     
%     
%     ub = 0.1;
% %     trace(PhiImmn*Xsdp) <= trace(PhiRemn*Xsdp)*tand(ub)
%     Phiub = PhiImmn - PhiRemn*tand(ub);
%     trace(Philb*Xsdp) <= 0;
    

    % Angle Reference
    
%     ref = -0.5;
%     
%     trace(PhiIm(:,:,4)*Xsdp) - trace(PhiRe(:,:,4)*Xsdp)*tand(ref) <= trace(PhiRe(:,:,4)*Xsdp)*gammasdp;
%     -trace(PhiIm(:,:,4)*Xsdp) + trace(PhiRe(:,:,4)*Xsdp)*tand(ref) <= trace(PhiRe(:,:,4)*Xsdp)*gammasdp;
    
%     Phimn = ea(:,6)*ea(:,4).';
%     PhiRemn = (0.5)*(Phimn + Phimn')
%     PhiImmn = (0.5/1j)*(Phimn - Phimn')
%     
%     trace(PhiImmn*Xsdp) <= gammasdp;
%     -trace(PhiImmn*Xsdp) <= gammasdp;

%     [PhiRemn j*PhiImmn;
%         j*PhiImmn PhiRemn]

%     [alphasdp*PhiRemn j*PhiImmn;
%         j*PhiImmn PhiRemn] >= 0

    % Angle Difference Minimization
    
%     Phimn = ea(:,4)*ea(:,6).';
%     PhiRemn = 0.5*(Phimn + Phimn');
%     PhiImmn = (0.5/1j)*(Phimn - Phimn');
    
%     [-gammasdp trace(PhiRemn*Xsdp);
%         trace(PhiRemn*Xsdp) 1] <= 0;

%     [gammasdp trace(PhiImmn*Xsdp);
%         trace(PhiImmn*Xsdp) 1] >= 0;

%     [gammasdp trace(PhiRemn*Xsdp) trace(PhiImmn*Xsdp);
%         trace(PhiRemn*Xsdp) -1 0;
%         trace(PhiImmn*Xsdp) 0 1] <= 0;

    
    % Feeder Head Voltage
    Xsdp(1,1) == V0(1)*V0(1)';
%     Xsdp >= 0;
    Xsdp == hermitian_semidefinite(nnode);
    
cvx_end;

%%

% [U, S, V] = svd(Xsdp);
% Xsdp = U(:,1)*S(1,1)*V(:,1)';

wsdp = zeros(3,nnode);
for k1 = 1:nnode
    
%     Vmagsdp(:,k1) = sqrt(diag(Xsdp(3*k1-2:3*k1,3*k1-2:3*k1)));


    %Phase a
    if wmaxpu(1,k1) == 0
        
    elseif wmaxpu(1,k1) > 0
        wsdp(1,k1) = trace(PhiPa(:,:,k1)*Xsdp) - real(spu(1,k1)) + 1j*(trace(PhiQa(:,:,k1)*Xsdp) - imag(spu(1,k1)) + caps.cappu(1,k1));
%         wsdp(1,k1) = trace(PhiPa(:,:,k1)*Xsdp) - real(spu(1,k1))*(aPQ(1,k1) + aZ(1,k1)*trace(PhiVa(:,:,k1)*Xsdp)) + ...
%         1j*(trace(PhiQa(:,:,k1)*Xsdp) - imag(spu(1,k1))*(aPQ(1,k1) + aZ(1,k1)*trace(PhiVa(:,:,k1)*Xsdp)) + caps.cappu(1,k1));
    end
    
%     %Phase b
%     if wmaxpu(2,k1) == 0
%         
%     elseif wmaxpu(2,k1) > 0
%         wsdp(2,k1) = trace(PhiPb(:,:,k1)*Xsdp) - real(spu(2,k1)) + 1j*(trace(PhiQb(:,:,k1)*Xsdp) - imag(spu(2,k1)) + caps.cappu(2,k1));
% %         wsdp(2,k1) = trace(PhiPb(:,:,k1)*Xsdp) - real(spu(2,k1))*(aPQ(2,k1) + aZ(2,k1)*trace(PhiVb(:,:,k1)*Xsdp)) + ...
% %         1j*(trace(PhiQb(:,:,k1)*Xsdp) - imag(spu(2,k1))*(aPQ(2,k1) + aZ(2,k1)*trace(PhiVb(:,:,k1)*Xsdp)) + caps.cappu(2,k1));
%     end
%     
%     %Phase c
%     if wmaxpu(3,k1) == 0
%         
%     elseif wmaxpu(3,k1) > 0
%         wsdp(3,k1) = trace(PhiPc(:,:,k1)*Xsdp) - real(spu(3,k1)) + 1j*(trace(PhiQc(:,:,k1)*Xsdp) - imag(spu(3,k1)) + caps.cappu(3,k1));
% %         wsdp(3,k1) = trace(PhiPc(:,:,k1)*Xsdp) - real(spu(3,k1))*(aPQ(3,k1) + aZ(3,k1)*trace(PhiVc(:,:,k1)*Xsdp)) + ...
% %         1j*(trace(PhiQc(:,:,k1)*Xsdp) - imag(spu(3,k1))*(aPQ(3,k1) + aZ(3,k1)*trace(PhiVc(:,:,k1)*Xsdp)) + caps.cappu(3,k1));
%     end
%     
%     Vangsdp(:,k1) = Xsdp(1,(3*k1)-2:3*k1)';
%     
%     Vsdp(:,k1) = Xsdp((3*k1)-2:3*k1,1);
    
end

% Vangsdp = 180/pi*angle(Vangsdp);

Vsdp = Xsdp(:,1);
Vmagsdp = abs(Xsdp);
Vangsdp = 180/pi*angle(Vsdp);


sdpsol.Xsdp = Xsdp;
% sdpsol.Ysdp = Ysdp;
sdpsol.Zsdp = Zsdp;
sdpsol.Vmagsdp = Vmagsdp;
sdpsol.Vangsdp = Vangsdp;
sdpsol.Vsdp = Vsdp;
sdpsol.wsdp = wsdp;
sdpsol.cvxstatus = cvx_status;