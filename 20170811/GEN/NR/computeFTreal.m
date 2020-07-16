function FT = computeFTreal(X,feeder,nodes,lines,configs,loads,caps,controllers,slackidx,slackVnom)

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

FTSUBV = [X(2*slackidx-1) - real(slackVnom(1));
    X(2*slackidx) - imag(slackVnom(1));
    X(2*nnode + 2*slackidx-1) - real(slackVnom(2));
    X(2*nnode + 2*slackidx) - imag(slackVnom(2));
    X(2*2*nnode + 2*slackidx-1) - real(slackVnom(3));
    X(2*2*nnode + 2*slackidx) - imag(slackVnom(3))];

FTKVL = zeros(2*3*nline,1);
for ph = 1:3
    for k1 = 1:nline
        
        idxre = 2*(ph-1)*nline + 2*k1-1;
        idxim = 2*(ph-1)*nline + 2*k1;
               
        idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1-1;
        idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1;

        if LPH(ph,k1) == 0

            FTKVL(idxre) = X(idxCmn);
            FTKVL(idxim) = X(idxDmn);

        elseif LPH(ph,k1) == 1
            
            idxAmTx = 2*(ph-1)*nnode + 2*TXnum(k1)-1;
            idxBmTx = 2*(ph-1)*nnode + 2*TXnum(k1);
            
            idxAmRx = 2*(ph-1)*nnode + 2*RXnum(k1)-1;            
            idxBmRx = 2*(ph-1)*nnode + 2*RXnum(k1);
            
            idxCmna = 2*3*nnode + 2*k1-1;
            idxDmna = 2*3*nnode + 2*k1;
            
            idxCmnb = 2*3*nnode + 2*nline + 2*k1-1;
            idxDmnb = 2*3*nnode + 2*nline + 2*k1;
            
            idxCmnc = 2*3*nnode + 2*2*nline + 2*k1-1;
            idxDmnc = 2*3*nnode + 2*2*nline + 2*k1;

            FTKVL(idxre) = X(idxAmTx) - X(idxAmRx) ...
                - real(FZpu(ph,1,k1))*X(idxCmna) + imag(FZpu(ph,1,k1))*X(idxDmna) ...
                - real(FZpu(ph,2,k1))*X(idxCmnb) + imag(FZpu(ph,2,k1))*X(idxDmnb) ...
                - real(FZpu(ph,3,k1))*X(idxCmnc) + imag(FZpu(ph,3,k1))*X(idxDmnc);
            if ph == 1 && k1 == 1
                
            end
            
            FTKVL(idxim) = X(idxBmTx) - X(idxBmRx) ...
                - real(FZpu(ph,1,k1))*X(idxDmna) - imag(FZpu(ph,1,k1))*X(idxCmna) ...
                - real(FZpu(ph,2,k1))*X(idxDmnb) - imag(FZpu(ph,2,k1))*X(idxCmnb) ...
                - real(FZpu(ph,3,k1))*X(idxDmnc) - imag(FZpu(ph,3,k1))*X(idxCmnc);
            
        end
        
    end
end

FTKCL = zeros(2*3*(nnode-1),1);
for ph = 1:3
    for k1 = 2:nnode
        
        idxre = 2*(ph-1)*nnode + 2*(k1-1)-1;
        idxim = 2*(ph-1)*nnode + 2*(k1-1);
        
        idxAm = 2*(ph-1)*nnode + 2*k1-1;
        idxBm = 2*(ph-1)*nnode + 2*k1;
        
        if NPH(ph,k1) == 0
            
            FTKCL(idxre) = X(idxAm);
            FTKCL(idxim) = X(idxBm);
            
        elseif NPH(ph,k1) == 1
            
            FTKCL(idxre) = 0;
            FTKCL(idxim) = 0;
            
            for k2 = 1:length(inmat(:,k1))
                
                if inmat(k2,k1) ~= 0
                    
                    idxClm = 2*3*nnode + 2*(ph-1)*nline + 2*inmat(k2,k1)-1;
                    idxDlm = 2*3*nnode + 2*(ph-1)*nline + 2*inmat(k2,k1);
                    
                    FTKCL(idxre) = FTKCL(idxre) + X(idxAm)*X(idxClm) + X(idxBm)*X(idxDlm);                    
                    FTKCL(idxim) = FTKCL(idxim) - X(idxAm)*X(idxDlm) + X(idxBm)*X(idxClm);
                    
                end
                
            end
            FTKCL(idxre) = FTKCL(idxre) - real(spu(ph,k1))*(APQ(ph,k1) + AZ(ph,k1)*(X(idxAm)^2 + X(idxBm)^2)) - real(wpu(ph,k1));
            FTKCL(idxim) = FTKCL(idxim) - imag(spu(ph,k1))*(APQ(ph,k1) + AZ(ph,k1)*(X(idxAm)^2 + X(idxBm)^2)) - imag(wpu(ph,k1)) + cappu(ph,k1);
            
            for k2 = 1:length(outmat(:,k1))
                
                if outmat(k2,k1) ~= 0
                    
                    idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*outmat(k2,k1)-1;
                    idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*outmat(k2,k1);
                    
                    FTKCL(idxre) = FTKCL(idxre) - X(idxAm)*X(idxCmn) - X(idxBm)*X(idxDmn);
                    FTKCL(idxim) = FTKCL(idxim) + X(idxAm)*X(idxDmn) - X(idxBm)*X(idxCmn);
                    
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