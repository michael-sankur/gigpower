function [VNRA1, STXNRA1, SRXNRA1] = NR3A1(network,Vopt,Iopt,subidx,Vnom,tol)

if ~exist('tol','var')
  tol=1e-9;
end

nnode = network.nodes.nnode;

nline = network.lines.nline;

nvar = nnode + nnode + nline + nline;

% XNR = []
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

XNRA1 = [1*ones(nnode,1);
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


FTA1 = 1e99;
iter = 0;
while max(abs(FTA1)) > tol
    
%     [a b] = max(abs(FTA1))
%     FTA1(b)

    FTA1 = computeFTrealA1(XNRA1,network,subidx,Vnom);

    JTA1 = computeJTrealA1(XNRA1,network,subidx);
    
    size(FTA1);
    size(JTA1);
    rank(JTA1);
    
%     eig(JTA1*JTA1.');
    
%     abs(eig(JTA1))
    
%     inv(JTA1);

%     eig(JTA1)

    JTA1;
    
    if size(JTA1,1) >= size(JTA1,2)
%         disp 'tall';
        XNRA1 = XNRA1 - inv(JTA1.'*JTA1)*JTA1.'*FTA1;
    elseif size(JTA1,1) < size(JTA1,2)
        rank(JTA1*JTA1.');
        XNRA1 = XNRA1 - JTA1.'*inv(JTA1*JTA1.')*FTA1;
    end
    
%     XNRA1 = XNRA1 - inv(JTA1)*FTA1;

    XNRA1;
           
    iter = iter + 1;
    
%     pause
        
end

FTA1;

XX(1,:) = XNRA1(1:nvar)';
XX(2,:) = XNRA1(nvar + 1:2*nvar)';
XX(3,:) = XNRA1(2*nvar + 1:end)';

ENRA1 = XX(:,1:nnode); XX(:,1:nnode) = [];
DNRA1 = XX(:,1:nnode); XX(:,1:nnode) = [];
PNRA1 = XX(:,1:nline); XX(:,1:nline) = [];
QNRA1 = XX(:,1:nline); clear XX;

VNRA1 = sqrt(ENRA1).*exp(j*DNRA1);

STXNRA1 = PNRA1 + 1j*QNRA1;
SRXNRA1 = PNRA1 + 1j*QNRA1;
