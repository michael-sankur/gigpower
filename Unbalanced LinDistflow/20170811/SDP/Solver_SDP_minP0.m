function [sdpsol,SDPMAT] = Solver_SDP_TX_minP0(feeder, nodes, lines, configs, loads, caps, controllers, sim)
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
    Hfbs1 = sim.Hfbs;
end

nvar = nnode + nnode + nline + nline + nnode + nnode; % number of variables per phase

%% Construct Y matrix

Y = zeros(3*nnode,3*nnode);

for k1 = 1:nline
    
    txn = lines.TXnum(k1);
    rxn = lines.RXnum(k1);
    
    Y(3*txn-2:3*txn,3*txn-2:3*txn) = Y(3*txn-2:3*txn,3*txn-2:3*txn) - FYpu(:,:,k1);
    Y(3*rxn-2:3*rxn,3*rxn-2:3*rxn) = Y(3*rxn-2:3*rxn,3*rxn-2:3*rxn) - FYpu(:,:,k1);
    Y(3*txn-2:3*txn,3*rxn-2:3*rxn) = FYpu(:,:,k1);
    Y(3*rxn-2:3*rxn,3*txn-2:3*txn) = FYpu(:,:,k1);
    
end

Y = Y;

% max(max(abs(Y - Y)))

%% Construct e vectors

e = eye(3*nnode);
for k1 = 1:nnode
    ea(:,k1) = e(:,3*k1-2);
    eb(:,k1) = e(:,3*k1-1);
    ec(:,k1) = e(:,3*k1);
end
% ea, eb, ec

%% Construct Phi and Y matrices

for k1 = 1:nnode
    
    PhiVa(:,:,k1) = ea(:,k1)*ea(:,k1).';
    PhiVb(:,:,k1) = eb(:,k1)*eb(:,k1).';
    PhiVc(:,:,k1) = ec(:,k1)*ec(:,k1).';
    
    Ya(:,:,k1) = Y'*PhiVa(:,:,k1);
    Yb(:,:,k1) = Y'*PhiVb(:,:,k1);
    Yc(:,:,k1) = Y'*PhiVc(:,:,k1);
    
%     Ya(:,:,k1) = PhiVa(:,:,k1)*Y;
%     Yb(:,:,k1) = PhiVb(:,:,k1)*Y;
%     Yc(:,:,k1) = PhiVc(:,:,k1)*Y;
    
    PhiIa(:,:,k1) = Y'*ea(:,k1)*ea(:,k1).'*Y;
    PhiIb(:,:,k1) = Y'*eb(:,k1)*eb(:,k1).'*Y;
    PhiIc(:,:,k1) = Y'*ec(:,k1)*ec(:,k1).'*Y;
       
    PhiPa(:,:,k1) = 0.5*(Ya(:,:,k1) + Ya(:,:,k1)');
    PhiPb(:,:,k1) = 0.5*(Yb(:,:,k1) + Yb(:,:,k1)');
    PhiPc(:,:,k1) = 0.5*(Yc(:,:,k1) + Yc(:,:,k1)');
    
%     PhiPa(:,:,k1) = PhiPa(:,:,k1) - aZ(1,k1)*real(spu(1,k1))*PhiVa(:,:,k1);
%     PhiPb(:,:,k1) = PhiPb(:,:,k1) - aZ(2,k1)*real(spu(2,k1))*PhiVb(:,:,k1);
%     PhiPc(:,:,k1) = PhiPc(:,:,k1) - aZ(3,k1)*real(spu(3,k1))*PhiVc(:,:,k1);
    
    PhiQa(:,:,k1) = -1j*0.5*(Ya(:,:,k1) - Ya(:,:,k1)');
    PhiQb(:,:,k1) = -1j*0.5*(Yb(:,:,k1) - Yb(:,:,k1)');
    PhiQc(:,:,k1) = -1j*0.5*(Yc(:,:,k1) - Yc(:,:,k1)');
    
%     PhiQa(:,:,k1) = PhiQa(:,:,k1) - aZ(1,k1)*imag(spu(1,k1))*PhiVa(:,:,k1);
%     PhiQb(:,:,k1) = PhiQb(:,:,k1) - aZ(2,k1)*imag(spu(2,k1))*PhiVb(:,:,k1);
%     PhiQc(:,:,k1) = PhiQc(:,:,k1) - aZ(3,k1)*imag(spu(3,k1))*PhiVc(:,:,k1);
    
end

%% Construct Phase Angle Matrices

% PhiVk0(:,:,1) = ea(:,1)*ea(:,tn)';
% PhiV0k(:,:,1) = ea(:,tn)*ea(:,1)';
% 
% PhiVk0(:,:,2) = eb(:,1)*eb(:,tn)';
% PhiV0k(:,:,2) = eb(:,tn)*eb(:,1)';
% 
% PhiVk0(:,:,3) = ec(:,1)*ec(:,tn)';
% PhiV0k(:,:,3) = ec(:,tn)*ec(:,1)';

% PhiAng = [];
% for ph = 1:3    
%     PhiAng(:,:,ph) = tand(delta_ref(ph) - delta0(ph))*(PhiVk0(:,:,ph) + PhiV0k(:,:,ph)) + j*(PhiVk0(:,:,ph) - PhiV0k(:,:,ph));
% end

% PhiAnga = [];
% PhiAngb = [];
% PhiAngc = [];
% for k1 = 2:nnode
%     if isnan(delta_ref(1,k1)) == 0
%         PhiAnga(:,:,k1) = tand(delta_ref(1,k1) - delta0(1))*(ea(:,1)*ea(:,k1)' + ea(:,k1)*ea(:,1)') + j*(ea(:,1)*ea(:,k1)' - ea(:,k1)*ea(:,1)');
%     end
%     if isnan(delta_ref(2,k1)) == 0
%         PhiAngb(:,:,k1) = tand(delta_ref(2,k1) - delta0(2))*(eb(:,1)*eb(:,k1)' + eb(:,k1)*eb(:,1)') + j*(eb(:,1)*eb(:,k1)' - eb(:,k1)*eb(:,1)');
%     end
%     if isnan(delta_ref(3,k1)) == 0
%         PhiAngc(:,:,k1) = tand(delta_ref(3,k1) - delta0(3))*(ec(:,1)*ec(:,k1)' + ec(:,k1)*ec(:,1)') + j*(ec(:,1)*ec(:,k1)' - ec(:,k1)*ec(:,1)');
%     end
% end

%%

SDPMAT.Y = Y;

SDPMAT.ea = ea;
SDPMAT.eb = eb;
SDPMAT.ec = ec;

SDPMAT.PhiVa = PhiVa;
SDPMAT.PhiVb = PhiVb;
SDPMAT.PhiVc = PhiVc;

SDPMAT.Ya = Ya;
SDPMAT.Yb = Yb;
SDPMAT.Yc = Yc;

SDPMAT.PhiPa = PhiPa;
SDPMAT.PhiPb = PhiPb;
SDPMAT.PhiPc = PhiPc;

SDPMAT.PhiQa = PhiQa;
SDPMAT.PhiQb = PhiQb;
SDPMAT.PhiQc = PhiQc;

% SDPMAT.PhiAnga = PhiAnga;
% SDPMAT.PhiAngb = PhiAngb;
% SDPMAT.PhiAngc = PhiAngc;

%%

wmaxpu(:,1) = 5*ones(3,1);

cvx_begin %sdp
%     cvx_solver sedumi;
    variable Xsdp(3*nnode,3*nnode) hermitian;
%     variable Ysdp(1,1) hermitian
    expression Zsdp;
    Zsdp = (trace(PhiVa(:,:,3)*Xsdp) + trace(PhiVb(:,:,3)*Xsdp) + trace(PhiVc(:,:,3)*Xsdp));
    maximize(Zsdp);
    subject to     
    
    for k1 = 1:nnode
        
        % Phase a
        if NPH(1,k1) == 0
            trace(PhiVa(:,:,k1)*Xsdp) == 0;
%             trace(PhiIa(:,:,k1)*Xsdp) == 0
%             trace(PhiPa(:,:,k1)*Xsdp) == 0
%             trace(PhiQa(:,:,k1)*Xsdp) == 0
        elseif NPH(1,k1) == 1
            
%             if wmaxpu(1,k1) == 0
%                 
%                 trace(PhiPa(:,:,k1)*Xsdp) == aPQ(1,k1)*real(spu(1,k1));
% 
%                 trace(PhiQa(:,:,k1)*Xsdp) == aPQ(1,k1)*imag(spu(1,k1)) - caps.cappu(1,k1);
%                 
%             elseif wmaxpu(1,k1) ~= 0
                
                trace(-PhiPa(:,:,k1)*Xsdp) >= -aPQ(1,k1)*real(spu(1,k1)) - wmaxpu(1,k1);
                trace(PhiPa(:,:,k1)*Xsdp) >= aPQ(1,k1)*real(spu(1,k1)) - wmaxpu(1,k1);

                trace(-PhiQa(:,:,k1)*Xsdp) >= -aPQ(1,k1)*imag(spu(1,k1)) + caps.cappu(1,k1) - wmaxpu(1,k1);
                trace(PhiQa(:,:,k1)*Xsdp) >= aPQ(1,k1)*imag(spu(1,k1)) - caps.cappu(1,k1) - wmaxpu(1,k1);
                
%             end

            trace(-PhiVa(:,:,k1)*Xsdp) >= -1.05^2;
            trace(PhiVa(:,:,k1)*Xsdp) >= 0.95^2;
            
        end

        % Phase b
        if NPH(2,k1) == 0
            trace(PhiVb(:,:,k1)*Xsdp) == 0;
%             trace(PhiIb(:,:,k1)*Xsdp) == 0
%             trace(PhiPb(:,:,k1)*Xsdp) == 0
%             trace(PhiQb(:,:,k1)*Xsdp) == 0
        elseif NPH(2,k1) == 1
            
%             if wmaxpu(2,k1) == 0
%                 
%                 trace(PhiPb(:,:,k1)*Xsdp) == aPQ(2,k1)*real(spu(2,k1));
% 
%                 trace(PhiQb(:,:,k1)*Xsdp) == aPQ(2,k1)*imag(spu(2,k1)) - caps.cappu(2,k1);
%                 
%             elseif wmaxpu(2,k1) ~= 0
            
                trace(-PhiPb(:,:,k1)*Xsdp) >= -aPQ(2,k1)*real(spu(2,k1)) - wmaxpu(2,k1);
                trace(PhiPb(:,:,k1)*Xsdp) >= aPQ(2,k1)*real(spu(2,k1)) - wmaxpu(2,k1);

                trace(-PhiQb(:,:,k1)*Xsdp) >= -aPQ(2,k1)*imag(spu(2,k1)) + caps.cappu(2,k1) - wmaxpu(2,k1);
                trace(PhiQb(:,:,k1)*Xsdp) >= aPQ(2,k1)*imag(spu(2,k1)) - caps.cappu(2,k1) - wmaxpu(2,k1);
            
%             end
            
            trace(-PhiVb(:,:,k1)*Xsdp) >= -1.05^2;
            trace(PhiVb(:,:,k1)*Xsdp) >= 0.95^2;
            
        end            

        % Phase c
        if NPH(3,k1) == 0
            trace(PhiVc(:,:,k1)*Xsdp) == 0;
%             trace(PhiIc(:,:,k1)*Xsdp) == 0
%             trace(PhiPc(:,:,k1)*Xsdp) == 0
%             trace(PhiQc(:,:,k1)*Xsdp) == 0
        elseif NPH(3,k1) == 1
            
%             if wmaxpu(3,k1) == 0
%                 
%                 trace(PhiPc(:,:,k1)*Xsdp) == aPQ(3,k1)*real(spu(3,k1));
% 
%                 trace(PhiQc(:,:,k1)*Xsdp) == aPQ(3,k1)*imag(spu(3,k1)) - caps.cappu(3,k1);
%                 
%             elseif wmaxpu(3,k1) ~= 0
            
                trace(-PhiPc(:,:,k1)*Xsdp) >= -aPQ(3,k1)*real(spu(3,k1)) - wmaxpu(3,k1);
                trace(PhiPc(:,:,k1)*Xsdp) >= aPQ(3,k1)*real(spu(3,k1)) - wmaxpu(3,k1);

                trace(-PhiQc(:,:,k1)*Xsdp) >= -aPQ(3,k1)*imag(spu(3,k1)) + caps.cappu(3,k1) - wmaxpu(3,k1);
                trace(PhiQc(:,:,k1)*Xsdp) >= aPQ(3,k1)*imag(spu(3,k1)) - caps.cappu(3,k1) - wmaxpu(3,k1);
                
%             end

            trace(-PhiVc(:,:,k1)*Xsdp) >= -1.05^2;
            trace(PhiVc(:,:,k1)*Xsdp) >= 0.95^2;
            
        end
        
    end
    
%     for k1 = 1:3*nnode
%         Xsdp(k1,k1) - Xsdp(k1,1)*Xsdp(1,k1) == 0
%     end
    
    % Feeder Head Voltage
    Xsdp(1,1) == V0(1)*V0(1)';
    Xsdp(1,2) == V0(1)*V0(2)';
    Xsdp(1,3) == V0(1)*V0(3)';
    Xsdp(2,1) == V0(2)*V0(1)';
    Xsdp(2,2) == V0(2)*V0(2)';
    Xsdp(2,3) == V0(2)*V0(3)';
    Xsdp(3,1) == V0(3)*V0(1)';
    Xsdp(3,2) == V0(3)*V0(2)';
    Xsdp(3,3) == V0(3)*V0(3)';
%     Xsdp >= 0;
    Xsdp == hermitian_semidefinite(3*nnode);
    
%     Ysdp == hermitian_semidefinite(1);
cvx_end;

%%

% [U, S, V] = svd(Xsdp);
% Xsdp = U(:,1)*S(1,1)*V(:,1)';

wsdp = zeros(3,nnode);
for k1 = 1:nnode
    
    Vmagsdp(:,k1) = sqrt(diag(Xsdp(3*k1-2:3*k1,3*k1-2:3*k1)));

    %Phase a
    if wmaxpu(1,k1) == 0
        
    elseif wmaxpu(1,k1) > 0
        wsdp(1,k1) = trace(PhiPa(:,:,k1)*Xsdp) - real(spu(1,k1)) + 1j*(trace(PhiQa(:,:,k1)*Xsdp) - imag(spu(1,k1)) + caps.cappu(1,k1));
%         wsdp(1,k1) = trace(PhiPa(:,:,k1)*Xsdp) - real(spu(1,k1))*(aPQ(1,k1) + aZ(1,k1)*trace(PhiVa(:,:,k1)*Xsdp)) + ...
%         1j*(trace(PhiQa(:,:,k1)*Xsdp) - imag(spu(1,k1))*(aPQ(1,k1) + aZ(1,k1)*trace(PhiVa(:,:,k1)*Xsdp)) + caps.cappu(1,k1));
    end
    
    %Phase b
    if wmaxpu(2,k1) == 0
        
    elseif wmaxpu(2,k1) > 0
        wsdp(2,k1) = trace(PhiPb(:,:,k1)*Xsdp) - real(spu(2,k1)) + 1j*(trace(PhiQb(:,:,k1)*Xsdp) - imag(spu(2,k1)) + caps.cappu(2,k1));
%         wsdp(2,k1) = trace(PhiPb(:,:,k1)*Xsdp) - real(spu(2,k1))*(aPQ(2,k1) + aZ(2,k1)*trace(PhiVb(:,:,k1)*Xsdp)) + ...
%         1j*(trace(PhiQb(:,:,k1)*Xsdp) - imag(spu(2,k1))*(aPQ(2,k1) + aZ(2,k1)*trace(PhiVb(:,:,k1)*Xsdp)) + caps.cappu(2,k1));
    end
    
    %Phase c
    if wmaxpu(3,k1) == 0
        
    elseif wmaxpu(3,k1) > 0
        wsdp(3,k1) = trace(PhiPc(:,:,k1)*Xsdp) - real(spu(3,k1)) + 1j*(trace(PhiQc(:,:,k1)*Xsdp) - imag(spu(3,k1)) + caps.cappu(3,k1));
%         wsdp(3,k1) = trace(PhiPc(:,:,k1)*Xsdp) - real(spu(3,k1))*(aPQ(3,k1) + aZ(3,k1)*trace(PhiVc(:,:,k1)*Xsdp)) + ...
%         1j*(trace(PhiQc(:,:,k1)*Xsdp) - imag(spu(3,k1))*(aPQ(3,k1) + aZ(3,k1)*trace(PhiVc(:,:,k1)*Xsdp)) + caps.cappu(3,k1));
    end
    
    Vangsdp(:,k1) = Xsdp(1,(3*k1)-2:3*k1)';
    
    Vsdp(:,k1) = Xsdp((3*k1)-2:3*k1,1);
    
end

Vangsdp = 180/pi*angle(Vangsdp);


sdpsol.Xsdp = Xsdp;
% sdpsol.Ysdp = Ysdp;
sdpsol.Zsdp = Zsdp;
sdpsol.Vmagsdp = Vmagsdp;
sdpsol.Vangsdp = Vangsdp;
sdpsol.Vsdp = Vsdp;
sdpsol.wsdp = wsdp;
sdpsol.cvxstatus = cvx_status;