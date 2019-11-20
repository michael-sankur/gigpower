function JT = computeJTrealA3(X,feeder,nodes,lines,configs,loads,caps,controllers,slackidx,slackVnom)

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

% controller parameters
wpu = controllers.wpu;

% Mmn and Nmn
Vfbs = [1*ones(1,nnode);
        1*exp(j*-120*pi/180)*ones(1,nnode);
        1*exp(j*120*pi/180)*ones(1,nnode)].*nodes.PH;
for k1 = 1:nline
    txnode = TXnum(k1);
    rxnode = RXnum(k1);
    gamman = Vfbs(:,txnode)*(1./Vfbs(:,txnode)).';
    if LPH(1,k1) == 0
        gamman(1,:) = 0; gamman(:,1) = 0;
    end
    if LPH(2,k1) == 0
        gamman(2,:) = 0; gamman(:,2) = 0;
    end
    if LPH(3,k1) == 0
        gamman(3,:) = 0; gamman(:,3) = 0;
    end
    gammaFZpu(:,:,k1) = gamman.*conj(FZpu(:,:,k1));
end
Mmn = real(gammaFZpu);
Nmn = imag(gammaFZpu);


nvar = 2*nnode + 2*nline;

% Substation Magnitude
% FTSUBMAG = [X(slackidx) - abs(slackVnom(1))^2;
%     X(nvar + slackidx) - abs(slackVnom(2))^2;
%     X(2*nvar + slackidx) - abs(slackVnom(3))^2];
JTSUBMAG = zeros(3,3*nvar);
JTSUBMAG(1,slackidx) = 1;
JTSUBMAG(2,nvar + slackidx) = 1;
JTSUBMAG(3,2*nvar + slackidx) = 1;

% Substation Angle
% FTSUBANG = [X(nnode + slackidx) - angle(slackVnom(1))^2;
%     X(nvar + nnode + slackidx) - angle(slackVnom(2))^2;
%     X(2*nvar + nnode + slackidx) - angle(slackVnom(3))^2];
JTSUBANG = zeros(3,3*nvar);
JTSUBANG(1,nnode + slackidx) = 1;
JTSUBANG(2,nvar + nnode + slackidx) = 1;
JTSUBANG(3,2*nvar + nnode + slackidx) = 1;

JTP0 = zeros(1,3*nvar);
JTQ0 = zeros(1,3*nvar);
JTMAG = zeros(1,3*nvar);
JTANG = zeros(1,3*nvar);
for ph = 1:3
    for k1 = 1:nline       
        
        if LPH(ph,k1) == 0
            
            idxPmn = (ph-1)*nvar + 2*nnode + k1;
            idxQmn = (ph-1)*nvar + 2*nnode + nline + k1;
            
%             FTP0(end+1,1) = X(idxPmn);
%             FTQ0(end+1,1) = X(idxQmn);

            JTP0(end+1,idxPmn) = 1;
            JTQ0(end+1,idxQmn) = 1;
            
        elseif LPH(ph,k1) == 1
            
            idxEm = (ph-1)*nvar + TXnum(k1);
            idxEn = (ph-1)*nvar + RXnum(k1);
            
            idxDm = (ph-1)*nvar + nnode + TXnum(k1);
            idxDn = (ph-1)*nvar + nnode + RXnum(k1);
            
            idxPmna = 2*nnode + k1;
            idxQmna = 2*nnode + nline + k1;
            
            idxPmnb = nvar + 2*nnode + k1;
            idxQmnb = nvar + 2*nnode + nline + k1;
            
            idxPmnc = 2*nvar + 2*nnode + k1;
            idxQmnc = 2*nvar + 2*nnode + nline + k1;
            
%             FTMAG(end+1,1) = -X(idxEm) + X(idxEn) ...
%                 + 2*(Mmn(ph,1,k1)*X(idxPmna) + Mmn(ph,2,k1)*X(idxPmnb) + Mmn(ph,3,k1)*X(idxPmnc)) ...
%                 - 2*(Nmn(ph,1,k1)*X(idxQmna) - Nmn(ph,2,k1)*X(idxQmnb) - Nmn(ph,3,k1)*X(idxQmnc));

            JTMAG(end+1,:) = zeros(1,3*nvar);
            JTMAG(end,idxEm) = -1;
            JTMAG(end,idxEn) = 1;
            JTMAG(end,idxPmna) = 2*Mmn(ph,1,k1);
            JTMAG(end,idxPmnb) = 2*Mmn(ph,2,k1);
            JTMAG(end,idxPmnc) = 2*Mmn(ph,3,k1);
            JTMAG(end,idxQmna) = -2*Nmn(ph,1,k1);
            JTMAG(end,idxQmnb) = -2*Nmn(ph,2,k1);
            JTMAG(end,idxQmnc) = -2*Nmn(ph,3,k1);
            
%             FTANG(end+1,1) = -X(idxDm) + X(idxDn) ...
%                 - Nmn(ph,1,k1)*X(idxPmna) - Nmn(ph,2,k1)*X(idxPmnb) - Nmn(ph,3,k1)*X(idxPmnc) ...
%                 - Mmn(ph,1,k1)*X(idxQmna) - Mmn(ph,2,k1)*X(idxQmnb) - Mmn(ph,3,k1)*X(idxQmnc); 

            JTANG(end+1,:) = zeros(1,3*nvar);
            JTANG(end,idxEm) = 0.5*sqrt(X(idxEn))*(-X(idxDm) + X(idxDn))/sqrt(X(idxEm));
            JTANG(end,idxEn) = 0.5*sqrt(X(idxEm))*(-X(idxDm) + X(idxDn))/sqrt(X(idxEn));
            JTANG(end,idxDm) = -sqrt(X(idxEm))*sqrt(X(idxEn));
            JTANG(end,idxDn) = sqrt(X(idxEm))*sqrt(X(idxEn));            
            JTANG(end,idxPmna) = -Nmn(ph,1,k1);
            JTANG(end,idxPmnb) = -Nmn(ph,2,k1);
            JTANG(end,idxPmnc) = -Nmn(ph,3,k1);
            JTANG(end,idxQmna) = -Mmn(ph,1,k1);
            JTANG(end,idxQmnb) = -Mmn(ph,2,k1);
            JTANG(end,idxQmnc) = -Mmn(ph,3,k1);
            
        end       
    end
end
JTP0 = JTP0(2:end,:);
JTQ0 = JTQ0(2:end,:);
JTMAG = JTMAG(2:end,:);
JTANG = JTANG(2:end,:);

% size(JTP0)
% size(JTQ0)
% size(JTMAG)
% size(JTANG)

JTE0 = zeros(1,3*nvar);
JTD0 = zeros(1,3*nvar);
JTPm = zeros(1,3*nvar);
JTQm = zeros(1,3*nvar);
for ph = 1:3
    for k1 = 2:nnode
        
        idxEm = (ph-1)*nvar + k1;
        
        idxDm = (ph-1)*nvar + nnode + k1;
        
        if NPH(ph,k1) == 0
                        
%             FTE0(end+1,1) = X(idxEm);
            
            JTE0(end+1,idxEm) = 1;
            
%             FTD0(end+1,1) = X(idxDm);

            JTD0(end+1,idxDm) = 1;
            
        elseif NPH(ph,k1) == 1
            
%             FTPm(end+1,1) = 0;
%             FTQm(end+1,1) = 0;
            
            JTPm(end+1,:) = zeros(1,3*nvar);
            JTQm(end+1,:) = zeros(1,3*nvar);
            
            for k2 = 1:length(inmat(:,k1))
                if inmat(k2,k1) ~= 0 && LPH(ph,inmat(k2,k1)) == 1
                                        
                    idxPlm = (ph-1)*nvar + 2*nnode + inmat(k2,k1);
                    idxQlm = (ph-1)*nvar + 2*nnode + nline + inmat(k2,k1);
                    
%                     FTPm(end) = FTPm(end) - X(idxPlm);
%                     FTQm(end) = FTQm(end) - X(idxQlm);
                    
                    JTPm(end,idxPlm) = -1;
                    JTQm(end,idxQlm) = -1;
                    
                end                
            end
            
%             FTPm(end) = FTPm(end) + real(spu(ph,k1))*(APQ(ph,k1) + AZ(ph,k1)*X(idxEm)) + real(wpu(ph,k1));
%             FTQm(end) = FTQm(end) + imag(spu(ph,k1))*(APQ(ph,k1) + AZ(ph,k1)*X(idxEm)) + imag(wpu(ph,k1)) - cappu(ph,k1);
            
            JTPm(end,idxEm) = real(spu(ph,k1))*AZ(ph,k1);
            JTQm(end,idxEm) = imag(spu(ph,k1))*AZ(ph,k1);
                        
            for k2 = 1:length(outmat(:,k1))
                if outmat(k2,k1) ~= 0 && LPH(ph,outmat(k2,k1)) == 1
                    
                    idxPmn = (ph-1)*nvar + 2*nnode + outmat(k2,k1);
                    idxQmn = (ph-1)*nvar + 2*nnode + nline + outmat(k2,k1);
                    
%                     FTPm(end) = JTPm(end) + X(idxPmn);
%                     FTQm(end) = JTQm(end) + X(idxQmn);
                    
                    JTPm(end,idxPmn) = 1;
                    JTQm(end,idxQmn) = 1;
                    
                end                
            end
            
        end
    end
end

JTE0 = JTE0(2:end,:);
JTD0 = JTD0(2:end,:);
JTPm = JTPm(2:end,:);
JTQm = JTQm(2:end,:);

% size(JTSUBMAG), rank(JTSUBMAG)
% size(JTSUBANG), rank(JTSUBANG)
% 
% size(JTP0), rank(JTP0)
% size(JTQ0), rank(JTQ0)
% size(JTMAG), rank(JTMAG)
% size(JTANG), rank(JTANG)
% 
% size(JTE0), rank(JTE0)
% size(JTD0), rank(JTD0)
% size(JTPm), rank(JTPm)
% size(JTQm), rank(JTQm)

% size([JTSUBMAG; JTE0; JTMAG; JTPm; JTQm; JTP0])
% rank([JTSUBMAG; JTE0; JTMAG; JTPm; JTQm; JTP0])
% 
% size([JTSUBANG; JTD0; JTANG;])
% rank([JTSUBANG; JTD0; JTANG;])
% 
% size([JTPm; JTQm; JTP0; JTQ0])
% rank([JTPm; JTQm; JTP0; JTQ0])

JT = [JTSUBMAG; JTSUBANG; JTP0; JTQ0; JTMAG; JTANG; JTE0; JTD0; JTPm; JTQm];

end