function [VNR, INR, STXNR, SRXNR] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,Vopt,Iopt,subidx,Vnom)

nnode = nodes.nnode;

nline = lines.nline;

XNR = [];
if isempty(Vopt) == 0
    for ph = 1:3    
        for k1 = 1:nnode
            XNR = [XNR;
                real(Vopt(ph,k1));
                imag(Vopt(ph,k1))];
        end
    end    
elseif isempty(Vopt) == 1
    for k1 = 1:nodes.nnode
        XNR = [XNR;
            real(Vnom(1));
            imag(Vnom(1));
            real(Vnom(2));
            imag(Vnom(2));
            real(Vnom(3));
            imag(Vnom(3))];
    end
end

if isempty(Iopt) == 0
    for ph = 1:3    
        for k1 = 1:nline
            XNR = [XNR;
                real(Iopt(ph,k1));
                imag(Iopt(ph,k1))];
        end
    end
elseif isempty(Iopt) == 1
    XNR = [XNR; 0.1*ones(6*lines.nline,1)];  
end


FT = 1e99;
iter = 0;
while max(abs(FT)) > 1e-9

    FT = computeFTreal(XNR,feeder,nodes,lines,configs,loads,caps,controllers,subidx,Vnom);

    JT = computeJTreal(XNR,feeder,nodes,lines,configs,loads,caps,controllers,subidx);
    
    if size(JT,1) >= size(JT,2)        
        XNR = XNR - inv(JT.'*JT)*JT.'*FT;
        ((eig(JT.'*JT)));
    end
       
    iter = iter + 1;
        
end

for k1 = 2:2:3*2*nnode
    VNR(k1/2) = XNR(k1-1) + 1j*XNR(k1);
end
VNR = [VNR(1:nnode);
    VNR(nnode+1:2*nnode);
    VNR(2*nnode+1:end)];

for k1 = 2:2:3*2*nline
    INR(k1/2) = XNR(3*2*nnode + k1-1) + 1j*XNR(3*2*nnode + k1);
end
INR = [INR(1:nline);
    INR(nline+1:2*nline);
    INR(2*nline+1:end)];

for k1 = 1:nline
%     if k1 >= 2
        STXNR(:,k1) = VNR(:,lines.TXnum(k1)).*conj(INR(:,k1));
%     end
    SRXNR(:,k1) = VNR(:,lines.RXnum(k1)).*conj(INR(:,k1));    
end
STXNR;
SRXNR;