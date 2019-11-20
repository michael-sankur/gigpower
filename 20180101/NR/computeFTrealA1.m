function FT = computeFTrealA1(X,network,slackidx,slackVnom)

% node parameters
nnode = network.nodes.nnode;
NPH = network.nodes.PH;
inmat = network.nodes.inmat;
outmat = network.nodes.outmat;

% line paramters
nline = network.lines.nline;
LPH = network.lines.PH;
TXnum = network.lines.TXnum;
RXnum = network.lines.RXnum;
FZpu = network.lines.FZpu;

% load parameters
spu = network.loads.spu;
APQ = network.loads.aPQ;
AI = network.loads.aI;
AZ = network.loads.aZ;

% capacitor paramters
cappu = network.caps.cappu;

% controller parameters
wpu = network.cons.wpu;

% Mmn and Nmn
Vfbs = [1*ones(1,nnode);
        1*exp(j*-120*pi/180)*ones(1,nnode);
        1*exp(j*120*pi/180)*ones(1,nnode)].*network.nodes.PH;
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


nvar = 2*nnode + nline + nline;

% Substation Magnitude
FTSUBMAG = [X(slackidx) - abs(slackVnom(1))^2;
    X(nvar + slackidx) - abs(slackVnom(2))^2;
    X(2*nvar + slackidx) - abs(slackVnom(3))^2];

% Substation Angle
FTSUBANG = [X(nnode + slackidx) - angle(slackVnom(1));
    X(nvar + nnode + slackidx) - angle(slackVnom(2));
    X(2*nvar + nnode + slackidx) - angle(slackVnom(3))];

FTP0 = [];
FTQ0 = [];
FTMAG = [];
FTANG = [];
for ph = 1:3
    for k1 = 1:nline
        
        if LPH(ph,k1) == 0
            
            idxPmn = (ph-1)*nvar + 2*nnode + k1;
            idxQmn = (ph-1)*nvar + 2*nnode + nline + k1;
            
            FTP0(end+1,1) = X(idxPmn);
            FTQ0(end+1,1) = X(idxQmn);
            
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
            
            FTMAG(end+1,1) = -X(idxEm) + X(idxEn) ...
                + 2*(Mmn(ph,1,k1)*X(idxPmna) + Mmn(ph,2,k1)*X(idxPmnb) + Mmn(ph,3,k1)*X(idxPmnc)) ...
                - 2*(Nmn(ph,1,k1)*X(idxQmna) + Nmn(ph,2,k1)*X(idxQmnb) + Nmn(ph,3,k1)*X(idxQmnc));
            
            FTANG(end+1,1) = -X(idxDm) + X(idxDn) ...
                - Nmn(ph,1,k1)*X(idxPmna) - Nmn(ph,2,k1)*X(idxPmnb) - Nmn(ph,3,k1)*X(idxPmnc) ...
                - Mmn(ph,1,k1)*X(idxQmna) - Mmn(ph,2,k1)*X(idxQmnb) - Mmn(ph,3,k1)*X(idxQmnc);            
            
        end       
    end
end


FTE0 = [];
FTD0 = [];
FTPm = [];
FTQm = [];
for ph = 1:3
    for k1 = 2:nnode
        
        idxEm = (ph-1)*nvar + k1;
        
        idxDm = (ph-1)*nvar + nnode + k1;
        
        if NPH(ph,k1) == 0
                        
            FTE0(end+1,1) = X(idxEm);
            
            FTD0(end+1,1) = X(idxDm);
            
        elseif NPH(ph,k1) == 1
            
            FTPm(end+1,1) = 0;
            FTQm(end+1,1) = 0;
            
            for k2 = 1:length(inmat(:,k1))                
                if inmat(k2,k1) ~= 0
                    
                    idxPlm = (ph-1)*nvar + 2*nnode + inmat(k2,k1);
                    idxQlm = (ph-1)*nvar + 2*nnode + nline + inmat(k2,k1);
                    
                    FTPm(end) = FTPm(end) - X(idxPlm);
                    FTQm(end) = FTQm(end) - X(idxQlm);
                    
                end                
            end
            
            FTPm(end) = FTPm(end) + real(spu(ph,k1))*(APQ(ph,k1) + AZ(ph,k1)*X(idxEm)) + real(wpu(ph,k1));
            FTQm(end) = FTQm(end) + imag(spu(ph,k1))*(APQ(ph,k1) + AZ(ph,k1)*X(idxEm)) + imag(wpu(ph,k1)) - cappu(ph,k1);
            
            for k2 = 1:length(outmat(:,k1))
                if outmat(k2,k1) ~= 0
                    
                    idxPmn = (ph-1)*nvar + 2*nnode + outmat(k2,k1);
                    idxQmn = (ph-1)*nvar + 2*nnode + nline + outmat(k2,k1);
                    
                    FTPm(end) = FTPm(end) + X(idxPmn);
                    FTQm(end) = FTQm(end) + X(idxQmn);
                    
                end                
            end           
        end        
    end
end

% FTSUBMAG
% FTSUBANG
% 
% FTP0
% FTQ0
% FTMAG
% FTANG
% 
% FTE0
% FTPm
% FTQm

% size(FTSUBMAG)
% size(FTSUBANG)
% 
% size(FTP0)
% size(FTQ0)
% size(FTMAG)
% size(FTANG)
% 
% size(FTE0)
% size(FTD0)
% size(FTPm)
% size(FTQm)


FT = [FTSUBMAG; FTSUBANG; FTP0; FTQ0; FTMAG; FTANG; FTE0; FTD0; FTPm; FTQm];

end