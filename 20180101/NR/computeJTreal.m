function JT = computeJTreal(XNR,network,slacknode)

% Michael Sankur - msankur@lbl.gov
% 2018.06.01

% This function computes the Jacobian of the residuals for the
% Newton-Raphson power flow solver.

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
% slacknode - index of slack node

% OUTPUT(S)
% JT - Partial derivatives of the residuals for power flow equations,
% composed of three parts - see near end of script
% JTSUBV - Partial derivatives of the residuals of slackbus real and
% imaginary voltage equation components
% JTKVL - Partial derivatives of the residuals of KVL real and imaginary
% equation components
% JTKCL - Partial derivatives of the residuals of KCL real and imaginary
% equation components

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

% node parameters
nnode = nodes.nnode;
NPH = nodes.PH;
inmat = nodes.inmat;
outmat = nodes.outmat;

% line paramters
nline = lines.nline;
LPH = lines.PH;
TXnum = lines.TXnum;
RXnum = lines.RXnum;
FZpu = lines.FZpu;

% load parameters
spu = loads.spu;
APQ = loads.aPQ;
AI = loads.aI;
AZ = loads.aZ;

% capacitor paramters
cappu = caps.cappu;

% Jacobian for slack node voltage
JSUBV = zeros(6,2*3*(nnode + nline));
for ph = 1:3
    
    % indexes of the real and imaginary components of the residual
    idxre = 2*ph-1;
    idxim = 2*ph;
    
    % indexes of the real and imaginary components voltage
    idxAm = 2*(ph-1)*nnode + 2*slacknode-1;
    idxBm = 2*(ph-1)*nnode + 2*slacknode;
    
    JSUBV(idxre,idxAm) = 1;
    JSUBV(idxim,idxBm) = 1;
    
end

% Jacobian for KVL across lines (m,n)
JKVL = zeros(2*3*nline,2*3*(nnode + nline));
for ph = 1:3
    for k1 = 1:nline
        
        % indexes of the real and imaginary components KVL residual
        idxre = 2*(ph-1)*nline + 2*k1-1;
        idxim = 2*(ph-1)*nline + 2*k1;
        
        % if phase does not exist on line
        % I_mn^phi = C_mn^phi + j D_mn^phi = 0
        if LPH(ph,k1) == 0
            
            % indexes of C_mn^phi and D_mn^phi
            idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1-1;
            idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1;

            % set derivatives of residuals for KVL
            JKVL(idxre,idxCmn) = 1;
            JKVL(idxim,idxDmn) = 1;

        % if phase does exist on line
        % V_m^phi = V_n^phi + sum_{psi} Z_{mn}^{phi psi} I_{mn}^psi
        % real: A_m^phi = A_n^phi + sum_{psi} r_{mn}^{phi psi} C_{mn}^psi - x_{mn}^{phi psi} D_{mn}^psi
        % imag: B_m^phi = B_n^phi + sum_{psi} r_{mn}^{phi psi} D_{mn}^psi + x_{mn}^{phi psi} C_{mn}^psi
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
            
            % set partial derivatives of residuals for KVL
            
            % real: A_m^phi - A_n^phi - sum_{psi} r_{mn}^{phi,psi} C_{mn}^psi - x_{mn}^{phi,psi} D_{mn}^psi = 0
            
            % derivatives of real KVl residual with respect to real component of node voltage
            JKVL(idxre,idxAmTx) = 1;
            JKVL(idxre,idxAnRx) = -1;
            
            % derivatives of real KVl residual with respect to real component of line current for phases a,b,c
            JKVL(idxre,idxCmna) = -real(FZpu(ph,1,k1));
            JKVL(idxre,idxCmnb) = -real(FZpu(ph,2,k1));
            JKVL(idxre,idxCmnc) = -real(FZpu(ph,3,k1));
            
            % derivatives of real KVl residual with respect to imag component of line current for phases a,b,c
            JKVL(idxre,idxDmna) = imag(FZpu(ph,1,k1));
            JKVL(idxre,idxDmnb) = imag(FZpu(ph,2,k1));
            JKVL(idxre,idxDmnc) = imag(FZpu(ph,3,k1));
            
            % imag: B_m^phi - B_n^phi - sum_{psi} r_{mn}^{phi,psi} D_{mn}^psi + x_{mn}^{phi,psi} C_{mn}^psi = 0
            
            % derivatives of real KVl residual with respect to imag component of node voltage
            JKVL(idxim,idxBmTx) = 1;
            JKVL(idxim,idxBnRx) = -1;
            
            % derivatives of imag KVl residual with respect to real component of line current for phases a,b,c
            JKVL(idxim,idxCmna) = -imag(FZpu(ph,1,k1));
            JKVL(idxim,idxCmnb) = -imag(FZpu(ph,2,k1));
            JKVL(idxim,idxCmnc) = -imag(FZpu(ph,3,k1));
            
            % derivatives of imag KVl residual with respect to imag component of line current for phases a,b,c
            JKVL(idxim,idxDmna) = -real(FZpu(ph,1,k1));
            JKVL(idxim,idxDmnb) = -real(FZpu(ph,2,k1));
            JKVL(idxim,idxDmnc) = -real(FZpu(ph,3,k1));
            
        end
                
    end
end

% Jacobian for KCL at node m
JKCL = zeros(2*3*(nnode-1),2*3*(nnode + nline));
for ph = 1:3
    for k1 = 2:nnode

        % indexes of real and imaginary KCL residual
        idxre = 2*(ph-1)*(nnode-1) + 2*(k1-1)-1;
        idxim = 2*(ph-1)*(nnode-1) + 2*(k1-1);

        % indexes of A_m^phi and B_m^phi
        idxAm = 2*(ph-1)*nnode + 2*k1-1;
        idxBm = 2*(ph-1)*nnode + 2*k1;

        % if phase does not exist at node, set V_m^phi = A_m^phi + j B_m^phi = 0
        if NPH(ph,k1) == 0

            JKCL(idxre,idxAm) = 1;
            JKCL(idxim,idxBm) = 1;

        % if phase does exist at node
        % sum_{l:(l,m) in Edges} V_m (I_lm^phi)^* = s_m^phi(V_m^phi) + w_m^phi - c_m^phi + sum_{n:(m,n) in Edges} V_m (I_mn^phi)^*
        elseif NPH(ph,k1) == 1

            % derivates of real KVL residual with respect to real and imag voltage components
            JKCL(idxre,idxAm) = -real(spu(ph,k1))*(AI(ph,k1)*XNR(idxAm)*(XNR(idxAm)^2 + XNR(idxBm)^2)^(-1/2) + AZ(ph,k1)*XNR(idxAm));            
            JKCL(idxre,idxBm) = -real(spu(ph,k1))*(AI(ph,k1)*XNR(idxBm)*(XNR(idxAm)^2 + XNR(idxBm)^2)^(-1/2) + AZ(ph,k1)*XNR(idxAm));

            % derivates of imag KVL residual with respect to real and imag voltage components
            JKCL(idxim,idxAm) = -imag(spu(ph,k1))*(AI(ph,k1)*XNR(idxAm)*(XNR(idxAm)^2 + XNR(idxBm)^2)^(-1/2) + AZ(ph,k1)*XNR(idxAm));
            JKCL(idxim,idxBm) = -imag(spu(ph,k1))*(AI(ph,k1)*XNR(idxBm)*(XNR(idxAm)^2 + XNR(idxBm)^2)^(-1/2) + AZ(ph,k1)*XNR(idxAm));

            % loop through incoming lines to node m - l:(l,m) in Edges
            for k2 = 1:length(inmat(:,k1))

                % incoming lines connected to node m
                if inmat(k2,k1) ~= 0

                    % indexes of I_lm^phi = C_mn^phi + j D_lm^phi
                    idxClm = 2*3*nnode + 2*(ph-1)*nline + 2*inmat(k2,k1)-1;
                    idxDlm = 2*3*nnode + 2*(ph-1)*nline + 2*inmat(k2,k1);

                    % derivaties of real KVL residual
                    JKCL(idxre,idxAm) = JKCL(idxre,idxAm) + XNR(idxClm);
                    JKCL(idxre,idxBm) = JKCL(idxre,idxBm) + XNR(idxDlm);                    
                    JKCL(idxre,idxClm) = XNR(idxAm);
                    JKCL(idxre,idxDlm) = XNR(idxBm);

                    % derivaties of imag KVL residual
                    JKCL(idxim,idxAm) = JKCL(idxim,idxAm) - XNR(idxDlm);                    
                    JKCL(idxim,idxBm) = JKCL(idxim,idxBm) + XNR(idxClm);                 
                    JKCL(idxim,idxClm) = XNR(idxBm);
                    JKCL(idxim,idxDlm) = -XNR(idxAm);

                end

            end

            % loop through outgoing lines from node m n:(m,n) in Edges
            for k2 = 1:length(outmat(:,k1))

                % outgoing lines connected to node m
                if outmat(k2,k1) ~= 0

                    % indexes of I_mn^phi = C_mn^phi + j D_mn^phi
                    idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*outmat(k2,k1)-1;
                    idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*outmat(k2,k1);

                    % derivaties of real KVL residual
                    JKCL(idxre,idxAm) = JKCL(idxre,idxAm) - XNR(idxCmn);
                    JKCL(idxre,idxBm) = JKCL(idxre,idxBm) - XNR(idxDmn);
                    JKCL(idxre,idxCmn) = -XNR(idxAm);
                    JKCL(idxre,idxDmn) = -XNR(idxBm);

                    % derivaties of imag KVL residual
                    JKCL(idxim,idxAm) = JKCL(idxim,idxAm) + XNR(idxDmn);
                    JKCL(idxim,idxBm) = JKCL(idxim,idxBm) - XNR(idxCmn);
                    JKCL(idxim,idxCmn) = -XNR(idxBm);
                    JKCL(idxim,idxDmn) = XNR(idxAm);

                end

            end

        end
    end
end

JT = [JSUBV; JKVL; JKCL];

end
