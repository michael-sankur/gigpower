function [VNRA2, STXNRA2, SRXNRA2] = NR3A3(feeder,nodes,lines,configs,loads,caps,controllers,Vopt,Iopt,subidx,Vnom)

nnode = nodes.nnode;

nline = lines.nline;

nvar = nnode + nnode + nline + nline;

% XNR = []h
% if isempty(Vopt) == 0
%     for ph
%     
%     
% elseif isempty(Vopt) == 1
%     
%     
% end

% XNR = [];
% if isempty(Vopt) == 0
%     for ph = 1:3
%         for k1 = 1:nnode
%             XNR = [XNR;
%                 real(Vopt(ph,k1));
%                 imag(Vopt(ph,k1))];
%         end
%     end    
% elseif isempty(Vopt) == 1
%     for k1 = 1:nodes.nnode
%         XNR = [XNR;
%             real(Vnom(1));
%             imag(Vnom(1));
%             real(Vnom(2));
%             imag(Vnom(2));
%             real(Vnom(3));
%             imag(Vnom(3))];
%     end
% end
% 
% if isempty(Iopt) == 0
%     for ph = 1:3    
%         for k1 = 1:nline
%             XNR = [XNR;
%                 real(Vopt(ph,k1));
%                 imag(Vopt(ph,k1))];
%         end
%     end
% elseif isempty(Iopt) == 1
%     XNR = [XNR; 0.1*ones(6*lines.nline,1)];  
% end

XNRA2 = [1*ones(nnode,1);
    0*ones(nnode,1);
    0.1*ones(nline,1);
    0.1*ones(nline,1);
    1*ones(nnode,1);
    -2*pi/3*ones(nnode,1);
    0.1*ones(nline,1);
    0.1*ones(nline,1);
    1*ones(nnode,1);
    2*pi/3*ones(nnode,1);
    0.1*ones(nline,1);
    0.1*ones(nline,1)];


FTA2 = 1e99;
iter = 0;
while max(abs(FTA2)) > 1e-9
    
%     [a b] = max(abs(FTA2))
%     FTA2(b)

    FTA2 = computeFTrealA2(XNRA2,feeder,nodes,lines,configs,loads,caps,controllers,subidx,Vnom);

    JTA2 = computeJTrealA2(XNRA2,feeder,nodes,lines,configs,loads,caps,controllers,subidx);
    
    size(FTA2);
    size(JTA2);
    rank(JTA2);
    
%     eig(JTA2*JTA2.');
    
%     abs(eig(JTA2))
    
%     inv(JTA2);

%     eig(JTA2)

    JTA2;
    
    if size(JTA2,1) >= size(JTA2,2)
%         disp 'tall';
        XNRA2 = XNRA2 - inv(JTA2.'*JTA2)*JTA2.'*FTA2;
    elseif size(JTA2,1) < size(JTA2,2)
        rank(JTA2*JTA2.');
        XNRA2 = XNRA2 - JTA2.'*inv(JTA2*JTA2.')*FTA2;
    end
    
%     XNRA2 = XNRA2 - inv(JTA2)*FTA2;

    XNRA2;
           
    iter = iter + 1;
    
%     pause
        
end

FTA2;

XX(1,:) = XNRA2(1:nvar)';
XX(2,:) = XNRA2(nvar + 1:2*nvar)';
XX(3,:) = XNRA2(2*nvar + 1:end)';

ENRA2 = XX(:,1:nnode); XX(:,1:nnode) = [];
DNRA2 = XX(:,1:nnode); XX(:,1:nnode) = [];
PNRA2 = XX(:,1:nline); XX(:,1:nline) = [];
QNRA2 = XX(:,1:nline); clear XX;

VNRA2 = sqrt(ENRA2).*exp(j*DNRA2);

STXNRA2 = PNRA2 + 1j*QNRA2;
SRXNRA2 = PNRA2 + 1j*QNRA2;
