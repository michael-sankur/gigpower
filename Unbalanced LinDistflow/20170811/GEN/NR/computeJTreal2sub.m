function JT = computeJTreal2sub(X,feeder,nodes,lines,configs,loads,caps,controllers,subidx)

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

JSUBV = zeros(6*length(subidx),2*3*(nnode + nline));
for k1 = 1:length(subidx)
    for ph = 1:3

        idxre = 6*(k1-1) + 2*ph-1;
        idxim = 6*(k1-1) + 2*ph;

        idxAm = 2*(ph-1)*nnode + 2*subidx(k1)-1;
        idxBm = 2*(ph-1)*nnode + 2*subidx(k1);

        JSUBV(idxre,idxAm) = 1;
        JSUBV(idxim,idxBm) = 1;

    end
end

JKVL = zeros(2*3*nline,2*3*(nnode + nline));
for ph = 1:3
    for k1 = 1:nline
        
        idxre = 2*(ph-1)*nline + 2*k1-1;
        idxim = 2*(ph-1)*nline + 2*k1;
               
        if LPH(ph,k1) == 0
            
            idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1-1;
            idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1;

            JKVL(idxre,idxCmn) = 1;
            JKVL(idxim,idxDmn) = 1;

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
            
            JKVL(idxre,idxAmTx) = 1;
            JKVL(idxre,idxAmRx) = -1;
            
            JKVL(idxre,idxCmna) = -real(FZpu(ph,1,k1));
            JKVL(idxre,idxCmnb) = -real(FZpu(ph,2,k1));
            JKVL(idxre,idxCmnc) = -real(FZpu(ph,3,k1));
            
            JKVL(idxre,idxDmna) = imag(FZpu(ph,1,k1));
            JKVL(idxre,idxDmnb) = imag(FZpu(ph,2,k1));
            JKVL(idxre,idxDmnc) = imag(FZpu(ph,3,k1));
            
            JKVL(idxim,idxBmTx) = 1;
            JKVL(idxim,idxBmRx) = -1;
            
            JKVL(idxim,idxCmna) = -imag(FZpu(ph,1,k1));
            JKVL(idxim,idxCmnb) = -imag(FZpu(ph,2,k1));
            JKVL(idxim,idxCmnc) = -imag(FZpu(ph,3,k1));
            
            JKVL(idxim,idxDmna) = -real(FZpu(ph,1,k1));
            JKVL(idxim,idxDmnb) = -real(FZpu(ph,2,k1));
            JKVL(idxim,idxDmnc) = -real(FZpu(ph,3,k1));
            
        end
                
    end
end

JKCL = zeros(2*3*(nnode-1),2*3*(nnode + nline));
for ph = 1:3
    for k1 = 2:nnode
        
        idxre = 2*(ph-1)*nnode + 2*(k1-1)-1;
        idxim = 2*(ph-1)*nnode + 2*(k1-1);
        
        idxAm = 2*(ph-1)*nnode + 2*k1-1;
        idxBm = 2*(ph-1)*nnode + 2*k1;
        
        if NPH(ph,k1) == 0
            
            JKCL(idxre,idxAm) = 1;
            JKCL(idxim,idxBm) = 1;
            
        elseif NPH(ph,k1) == 1
            
            JKCL(idxre,idxAm) = -2*real(spu(ph,k1))*AZ(ph,k1)*X(idxAm);            
            JKCL(idxre,idxBm) = -2*real(spu(ph,k1))*AZ(ph,k1)*X(idxBm);
            
            JKCL(idxim,idxAm) = -2*imag(spu(ph,k1))*AZ(ph,k1)*X(idxAm);
            JKCL(idxim,idxBm) = -2*imag(spu(ph,k1))*AZ(ph,k1)*X(idxBm);
            
            for k2 = 1:length(inmat(:,k1))
                
                if inmat(k2,k1) ~= 0
                    
                    idxClm = 2*3*nnode + 2*(ph-1)*nline + 2*inmat(k2,k1)-1;
                    idxDlm = 2*3*nnode + 2*(ph-1)*nline + 2*inmat(k2,k1);
                                        
                    JKCL(idxre,idxAm) = JKCL(idxre,idxAm) + X(idxClm);
                    JKCL(idxre,idxBm) = JKCL(idxre,idxBm) + X(idxDlm);                    
                    JKCL(idxre,idxClm) = X(idxAm);
                    JKCL(idxre,idxDlm) = X(idxBm);
                    
                    
                    JKCL(idxim,idxAm) = JKCL(idxim,idxAm) - X(idxDlm);                    
                    JKCL(idxim,idxBm) = JKCL(idxim,idxBm) + X(idxClm);                 
                    JKCL(idxim,idxClm) = X(idxBm);
                    JKCL(idxim,idxDlm) = -X(idxAm);
                    
                end
                
            end
                        
            for k2 = 1:length(outmat(:,k1))
                
                if outmat(k2,k1) ~= 0
                    
                    idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*outmat(k2,k1)-1;
                    idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*outmat(k2,k1);
                    
                    JKCL(idxre,idxAm) = JKCL(idxre,idxAm) - X(idxCmn);
                    JKCL(idxre,idxBm) = JKCL(idxre,idxBm) - X(idxDmn);
                    JKCL(idxre,idxCmn) = -X(idxAm);
                    JKCL(idxre,idxDmn) = -X(idxBm);
                    
                    JKCL(idxim,idxAm) = JKCL(idxim,idxAm) + X(idxDmn);
                    JKCL(idxim,idxBm) = JKCL(idxim,idxBm) - X(idxCmn);
                    JKCL(idxim,idxCmn) = -X(idxBm);
                    JKCL(idxim,idxDmn) = X(idxAm);
                    
                end
                
            end
            
        end
    end
end

JT = [JSUBV; JKVL; JKCL];
