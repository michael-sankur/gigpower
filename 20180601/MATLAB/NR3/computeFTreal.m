function FT = computeFTreal(XNR,network,slacknode,Vslack)

% Michael Sankur - msankur@lbl.gov
% 2018.06.01

% This function computes the residuals of unbalanced power flow equations
% (KVL and KCL)

% INPUT(S)
% XNR - current iteration of NR algorithm variable
% network - struct containing all pertinent the network information,
% including all other structs
% base - struct containing base values
% nodes - struct containing node parameters
% lines - struct containing line parameters
% loads - struct containing load parameters
% caps - struct containing capacitor parameters
% cons - struct containing controller parameters
% vvc - struct containing vvc parameters
% slacknode - index of slack node
% Vslack - voltage reference for slack node

% OUTPUT(S)
% FT - Residuals for power flow equations, composed of the following three
% parts - see end of function
% FTSUBV - residuals of slackbus real and imaginary voltage equation
% components
% FTKVL - residuals of KVL real and imaginary equation components
% FTKCL - residuals of KCL real and imaginary equation components

% slackidx is the node index of the slack bus, which is assigned a fixed
% voltage reference of slackVnom.

% Voltage and current are separated into their real and imaginary parts
% V_n^phi = A_n^phi + j B_n^phi
% I_n^phi = C_n^phi + j D_n^phi

% Voltage and current vectors for a single phase
% V^phi = [A_1^phi, B_1^phi, A_2^phi, B_2^phi, ... , A_n^phi, B_n^phi]
% I^phi = [C_1^phi, D_1^phi, C_2^phi, D_2^phi, ... , C_n^phi, D_n^phi]

% The NR algorithm variable
% X = [V^a V^b V^c I^a I^b I^c]

base = network.base;
nodes = network.nodes;
lines = network.lines;
configs = network.configs;
loads = network.loads;
caps = network.caps;
cons = network.cons;
vvc = network.vvc;

% node parameters
nnode = nodes.nnode;
NPH = nodes.PH;
nnode = nodes.nnode;
NPH = nodes.PH;
inlines = nodes.inlines;
innodes = nodes.innodes;
outlines = nodes.outlines;
outnodes = nodes.outnodes;

% line paramters
nline = lines.nline;
LPH = lines.PH;
TXnum = lines.TXnum;
RXnum = lines.RXnum;
FZpu = lines.FZpu;
FRpu = lines.FRpu;
FXpu = lines.FXpu;

% load parameters
spu = loads.spu;
APQ = loads.aPQ;
AI = loads.aI;
AZ = loads.aZ;

% capacitor paramters
cappu = caps.cappu;

% controller parameters
wpu = cons.wpu;

% vvc parameters
vvcpu = vvc.vvcpu;

% Residuals for slack node voltage
FTSUBV = [XNR(2*slacknode-1) - real(Vslack(1));
    XNR(2*slacknode) - imag(Vslack(1));
    XNR(2*nnode + 2*slacknode-1) - real(Vslack(2));
    XNR(2*nnode + 2*slacknode) - imag(Vslack(2));
    XNR(2*2*nnode + 2*slacknode-1) - real(Vslack(3));
    XNR(2*2*nnode + 2*slacknode) - imag(Vslack(3))];



% Residuals for KVL across line (m,n)
FTKVL = zeros(2*3*nline,1);
for ph = 1:3
    for k1 = 1:nline
        
        % indexes of real and imag parts of KVL equation for phase phi on line (m,n)
        idxre = 2*(ph-1)*nline + 2*k1-1;
        idxim = 2*(ph-1)*nline + 2*k1;
        
        % if phase does not exist on line
        % I_mn^phi = C_mn^phi + j D_mn^phi = 0
        if LPH(ph,k1) == 0
            
            % indexes of C_mn^phi and D_mn^phi
            idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1-1;
            idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1;

            % set residuals for KVL
            FTKVL(idxre) = XNR(idxCmn);
            FTKVL(idxim) = XNR(idxDmn);

        % if phase does exist on line
        % V_m^phi = V_n^phi + sum_{psi} Z_{mn}^{phi psi} I_{mn}^psi
        % Real: A_m^phi = A_n^phi + sum_{psi} r_{mn}^{phi psi} C_{mn}^psi - x_{mn}^{phi psi} D_{mn}^psi
        % Imag: B_m^phi = B_n^phi + sum_{psi} r_{mn}^{phi psi} D_{mn}^psi + x_{mn}^{phi psi} C_{mn}^psi
        elseif LPH(ph,k1) == 1
            
            % indexes of V_m^phi = A_m^phi +j B_m^phi for TX node
            idxAmTx = 2*(ph-1)*nnode + 2*TXnum(k1)-1;
            idxBmTx = 2*(ph-1)*nnode + 2*TXnum(k1);
            
            % indexes of V_n^phi = A_n^phi +j B_n^phi for RX node
            idxAnRx = 2*(ph-1)*nnode + 2*RXnum(k1)-1;            
            idxBnRx = 2*(ph-1)*nnode + 2*RXnum(k1);
            
            % indexes of I_mn^a = C_mn^a + j D_mn^a
            idxCmna = 2*3*nnode + 2*k1-1;
            idxDmna = 2*3*nnode + 2*k1;
            
            % indexes of I_mn^b = C_mn^b + j D_mn^b
            idxCmnb = 2*3*nnode + 2*nline + 2*k1-1;
            idxDmnb = 2*3*nnode + 2*nline + 2*k1;
            
            % indexes of I_mn^c = C_mn^c + j D_mn^c
            idxCmnc = 2*3*nnode + 2*2*nline + 2*k1-1;
            idxDmnc = 2*3*nnode + 2*2*nline + 2*k1;
            
            % set residuals for KVL
            
            % real: X_m^phi = A_m^phi - A_n^phi - sum_{psi} r_{mn}^{phi,psi} C_{mn}^psi - x_{mn}^{phi,psi} D_{mn}^psi
            FTKVL(idxre) = XNR(idxAmTx) - XNR(idxAnRx) ...
                - FRpu(ph,1,k1)*XNR(idxCmna) + FXpu(ph,1,k1)*XNR(idxDmna) ...
                - FRpu(ph,2,k1)*XNR(idxCmnb) + FXpu(ph,2,k1)*XNR(idxDmnb) ...
                - FRpu(ph,3,k1)*XNR(idxCmnc) + FXpu(ph,3,k1)*XNR(idxDmnc);
            
            % imag: Y_m^phi = B_m^phi - B_n^phi - sum_{psi} r_{mn}^{phi,psi} D_{mn}^psi + x_{mn}^{phi,psi} C_{mn}^psi
            FTKVL(idxim) = XNR(idxBmTx) - XNR(idxBnRx) ...
                - FRpu(ph,1,k1)*XNR(idxDmna) - FXpu(ph,1,k1)*XNR(idxCmna) ...
                - FRpu(ph,2,k1)*XNR(idxDmnb) - FXpu(ph,2,k1)*XNR(idxCmnb) ...
                - FRpu(ph,3,k1)*XNR(idxDmnc) - FXpu(ph,3,k1)*XNR(idxCmnc);
            
        end
        
    end
end

% Residuals for KCL at node m
% This algorithm assumes that the slack bus has a fixed voltage reference,
% and its power is "floating" and will be resolved. The slack bus is
% assumed to be the first node, which respresents the transmission line, or
% substation if the network configuration is as such - see note below
FTKCL = zeros(2*3*(nnode-1),1);
for ph = 1:3
    for k1 = 2:nnode
        
        % indexes of real and imag parts of KCL equation for node m
        idxre = 2*(ph-1)*(nnode-1) + 2*(k1-1)-1;
        idxim = 2*(ph-1)*(nnode-1) + 2*(k1-1);

        % indexes of A_m^phi and B_m^phi
        idxAm = 2*(ph-1)*nnode + 2*k1-1;
        idxBm = 2*(ph-1)*nnode + 2*k1;

        % if phase does not exist at node, set V_m^phi = A_m^phi + j B_m^phi = 0
        if NPH(ph,k1) == 0

            FTKCL(idxre) = XNR(idxAm);
            FTKCL(idxim) = XNR(idxBm);

        % if phase does exist at node
        % sum_{l:(l,m) in Edges} V_m (I_lm^phi)^* = s_m^phi(V_m^phi) + w_m^phi - c_m^phi + sum_{n:(m,n) in Edges} V_m (I_mn^phi)^*
        elseif NPH(ph,k1) == 1

            % initialize residual as zero
            FTKCL(idxre) = 0;
            FTKCL(idxim) = 0;

            % loop through incoming lines to node m - l:(l,m) in Edges
            for k2 = 1:length(inlines(:,k1))
                
                % incoming lines connected to node m
                if inlines(k2,k1) ~= 0

                    % indexes of I_lm^phi = C_lm^phi + j D_lm^phi
                    idxClm = 2*3*nnode + 2*(ph-1)*nline + 2*inlines(k2,k1)-1;
                    idxDlm = 2*3*nnode + 2*(ph-1)*nline + 2*inlines(k2,k1);

                    % sum_{l:(l,m) in Edges} A_m (I_lm^phi)^*
                    % real: A_m^phi C_lm^phi + B_m^phi D_lm^phi
                    % imag: -A_m^phi D_lm^phi + B_m^phi C_lm^phi
                    FTKCL(idxre) = FTKCL(idxre) + XNR(idxAm)*XNR(idxClm) + XNR(idxBm)*XNR(idxDlm);                    
                    FTKCL(idxim) = FTKCL(idxim) - XNR(idxAm)*XNR(idxDlm) + XNR(idxBm)*XNR(idxClm);

                end

            end

            % s_m^phi(V_m^phi) + w_m^phi - 1j*c_m^phi + 1j*vvc_m^phi
            % real: p_m^phi (A_{PQ,m}^phi + A_{I,m}^phi ((A_m^phi)^2 + (B_m^phi)^2)^(1/2) + A_{Z,m}^phi ((A_m^phi)^2 + (B_m^phi)^2)) - u_m^phi
            % imag: q_m^phi (A_{PQ,m}^phi + A_{I,m}^phi ((A_m^phi)^2 + (B_m^phi)^2)^(1/2) + A_{Z,m}^phi ((A_m^phi)^2 + (B_m^phi)^2)) + c_m^phi - v_m^phi - vvc_m^phi 
            FTKCL(idxre) = FTKCL(idxre) - real(spu(ph,k1))*(APQ(ph,k1) ...
                + AI(ph,k1)*(XNR(idxAm)^2 + XNR(idxBm)^2)^(1/2) ...
                + AZ(ph,k1)*(XNR(idxAm)^2 + XNR(idxBm)^2)) ...
                - real(wpu(ph,k1));
            FTKCL(idxim) = FTKCL(idxim) - imag(spu(ph,k1))*(APQ(ph,k1) ...
                + AI(ph,k1)*(XNR(idxAm)^2 + XNR(idxBm)^2)^(1/2) ...
                + AZ(ph,k1)*(XNR(idxAm)^2 + XNR(idxBm)^2)) ...
                + cappu(ph,k1) - imag(wpu(ph,k1)) - vvcpu(ph,k1);

            % loop through outgoing lines from node m - n:(m,n) in Edges
            for k2 = 1:length(outlines(:,k1))
                
                % outgoing lines connected to node m
                if outlines(k2,k1) ~= 0

                    % indexes of I_mn^phi = C_mn^phi + j D_mn^phi
                    idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*outlines(k2,k1)-1;
                    idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*outlines(k2,k1);

                    % sum_{n:(m,n) in Edges} A_m (I_mn^phi)^*
                    % real: -A_m^phi C_mn^phi - B_m^phi D_mn^phi
                    % imag: A_m^phi D_mn^phi - B_m^phi C_nm^phi
                    FTKCL(idxre) = FTKCL(idxre) - XNR(idxAm)*XNR(idxCmn) - XNR(idxBm)*XNR(idxDmn);
                    FTKCL(idxim) = FTKCL(idxim) + XNR(idxAm)*XNR(idxDmn) - XNR(idxBm)*XNR(idxCmn);

                end

            end

        end
    end
end

FTSUBV;
FTKVL;
FTKCL;

FT = [FTSUBV; FTKVL; FTKCL];

end

% Future versions may have ability to specify a node with fixed voltage,
% and a node with "floating" power which will be resolved by the algorithm.
% These nodes are may be the same or different.