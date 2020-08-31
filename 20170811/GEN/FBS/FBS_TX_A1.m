function [FBS] = FBS_TX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vref)

% Linear approximation
% No losses
% No Hmn
% Gamma based on nominal voltage
% Small angle approximation
% Assume voltage magnitudes = 1 for angle equation

% Michael Sankur - msankur@berkeley.edu
% 2016.06.03
% Adapted from function originally by Dan Arnold

% Node paramaters
nnode = nodes.nnode;

% Line parameters
nline = lines.nline;
FZpu = lines.FZpu;

% Load parameters
spu = loads.spu;
aPQ = loads.aPQ;
aI = loads.aI;
aZ = loads.aZ;

% Capacitor parameters
cappu = caps.cappu;

% Controller parameters
wpu = controllers.wpu;

% Find root, passthrough, junction and tail nodes
rn = 1;
pn = [];
jn = [1 2];
tn = [];
for k1 = 1:nodes.nnode
    if sum(nodes.inmat(:,k1) ~= 0) == 1 && sum(nodes.outmat(:,k1) ~= 0) == 0
        tn = [tn k1];
    end
    if sum(nodes.inmat(:,k1) ~= 0) == 1 && sum(nodes.outmat(:,k1) ~= 0) == 1
        pn = [pn k1];
    end
    if sum(nodes.inmat(:,k1) ~= 0) == 1 && sum(nodes.outmat(:,k1) ~= 0) >= 2
        jn = [jn k1];
    end
    
end
% rn, pn, jn, tn

% Initialize iteration
iter = 0;

gamma = Vref*(1./Vref).';
for k1 = 1:nline
    gammaFZpu(:,:,k1) = gamma.*conj(FZpu(:,:,k1));
end
Mmn = real(gammaFZpu);
Nmn = -imag(gammaFZpu);

E = zeros(3,nnode);
D = zeros(3,nnode);
P = zeros(3,nline);
Q = zeros(3,nline);

E0 = [1 1 1]';
Etest = [0 0 0]';

D0 = [0 240 120]';
Dtest = [0 0 0]';

Lmn = zeros(3,nline);
Hmn = zeros(3,nline);
V = ones(3,nnode);

while max(abs(Etest - E0)) >= 1e-6 && max(abs(Dtest - D0)) >= 1e-6
    
    E(:,1:2) = [E0 E0];
    D(:,1:2) = [D0 D0];
    % sweep forward, calculate voltages
    for k1 = 2:1:nnode
        % find children of current node
        lineidx = nodes.outmat(nodes.outmat(:,k1) ~= 0,k1);
        for k2 = 1:1:length(lineidx)
            kline = lineidx(k2);
            knode = lines.RXnum(kline);
%             kline = nodes.outmat(k2,k1);
            % calculate voltage at child node
            E(:,knode) = E(:,k1) - 2*(Mmn(:,:,kline)*P(:,kline) + Nmn(:,:,kline)*Q(:,kline));
            D(:,knode) = D(:,k1) - 180/pi*(Nmn(:,:,kline)*P(:,kline) - Mmn(:,:,kline)*Q(:,kline));
            E(nodes.PH == 0) = 0;
        end
        clear lineidx;
    end
    
    % recalculate voltage dependent load
    s = spu.*(aPQ + aI.*sqrt(E) + aZ.*E) - 1j*cappu + wpu;
    for k2 = nnode:-1:2
        inline = nodes.inmat(k2);
        outlines = nodes.outmat(nodes.outmat(:,k2) ~= 0,k2);
        P(:,inline) = real(s(:,k2)) + sum(P(:,outlines),2);
        Q(:,inline) = imag(s(:,k2)) + sum(Q(:,outlines),2);
    end
    P(lines.PH == 0) = 0; Q(lines.PH == 0) = 0;
    
    % TERMINAL NODES
%     disp 'Terminal Nodes'
    % loop through terminal nodes
    for k1 = length(tn):-1:1
%         disp('NEW TERMINAL NODE')
        cn = tn(k1);
        kline = nodes.inmat(cn);
        
        % recalculate voltage dependent load
        s = spu.*(aPQ + aI.*sqrt(E) + aZ.*E) - 1j*cappu + wpu;
        for k2 = nnode:-1:2
            inline = nodes.inmat(k2);
            outlines = nodes.outmat(nodes.outmat(:,k2) ~= 0,k2);
            P(:,inline) = real(s(:,k2)) + sum(P(:,outlines),2);
            Q(:,inline) = imag(s(:,k2)) + sum(Q(:,outlines),2);
        end
        P(lines.PH == 0) = 0; Q(lines.PH == 0) = 0;
        % calculate segment current
                
        % move up the feeder until a junction node is hit
        % assume each node only has 1 parent
        flag = true;
%         idx = find(FM(cn,:) == -1);
        idx = lines.TXnum(kline);
        while(flag)
            % calculate voltage at parent node only for phases existing at
            % current node. doing so for all phases will set a phase to
            % an erroneous value if it exists at the parent and not at the
            % current node
            for ph = 1:3
                if nodes.PH(ph,cn) == 1
                   E(ph,idx) = E(ph,cn) + 2*(Mmn(ph,:,kline)*P(:,kline) + Nmn(ph,:,kline)*Q(:,kline));
                   D(ph,idx) = D(ph,cn) + 180/pi*(Nmn(ph,:,kline)*P(:,kline) - Mmn(ph,:,kline)*Q(:,kline));
                   E(nodes.PH == 0) = 0;
                end
            end
            
            % recalculate voltage dependent load
            s = spu.*(aPQ + aI.*sqrt(E) + aZ.*E) - 1j*cappu + wpu;
            for k2 = nnode:-1:2
                inline = nodes.inmat(k2);
                outlines = nodes.outmat(nodes.outmat(:,k2) ~= 0,k2);
                P(:,inline) = real(s(:,k2)) + sum(P(:,outlines),2);
                Q(:,inline) = imag(s(:,k2)) + sum(Q(:,outlines),2);
            end
            P(lines.PH == 0) = 0; Q(lines.PH == 0) = 0;
            
            % check if idx is a junction node
            if(isempty(find(jn == idx)))
                % then this is not a junction node, move up path
%                 disp 'move on';
                cn = idx;
                kline = nodes.inmat(cn);
%                 idx = find(FM(cn,:) == -1);
                idx = lines.TXnum(nodes.inmat(cn));
            else
                % this is a junction node and we stop this branch
%                 disp 'end path, junction reached';
                flag = false;
            end
%             pause
        end
        % completed branches from terminal nodes to junction nodes
        % junction nodes are now terminal branches
    end
                
    % recalculate voltage dependent load
    s = spu.*(aPQ + aI.*sqrt(E) + aZ.*E) - 1j*cappu + wpu;
    for k2 = nnode:-1:2
        inline = nodes.inmat(k2);
        outlines = nodes.outmat(nodes.outmat(:,k2) ~= 0,k2);
        P(:,inline) = real(s(:,k2)) + sum(P(:,outlines),2);
        Q(:,inline) = imag(s(:,k2)) + sum(Q(:,outlines),2);
    end
    P(lines.PH == 0) = 0; Q(lines.PH == 0) = 0;
        
    % JUNCTION NODES
    % start from bottom-most junction node and repeat until a junction
    % node is hit, stop if junction node is 1
    % loop through junction nodes
%     disp 'Junction Nodes'
    for k1 = length(jn):-1:2
%         disp('NEW JUNCTION NODE')
        cn = jn(k1);
        kline = nodes.inmat(cn);

        % recalculate voltage dependent load
        s = spu.*(aPQ + aI.*sqrt(E) + aZ.*E) - 1j*cappu + wpu;
        for k2 = nnode:-1:2
            inline = nodes.inmat(k2);
            outlines = nodes.outmat(nodes.outmat(:,k2) ~= 0,k2);
            P(:,inline) = real(s(:,k2)) + sum(P(:,outlines),2);
            Q(:,inline) = imag(s(:,k2)) + sum(Q(:,outlines),2);
        end
        P(lines.PH == 0) = 0; Q(lines.PH == 0) = 0;
        
        % move up the feeder until a junction node is hit
        flag = true;
%         idx = find(FM(cn,:) == -1);
        idx = lines.TXnum(kline);
        while(flag)
            % calculate voltage at parent node only for phases existing at
            % current node. doing so for all phases will set a phase to
            % an erroneous value if it exists at the parent and not at the
            % current node
            for ph = 1:3
                if nodes.PH(ph,cn) == 1
                   E(ph,idx) = E(ph,cn) + 2*(Mmn(ph,:,kline)*P(:,kline) + Nmn(ph,:,kline)*Q(:,kline));
                   D(ph,idx) = D(ph,cn) + 180/pi*(Nmn(ph,:,kline)*P(:,kline) - Mmn(ph,:,kline)*Q(:,kline));
                   E(nodes.PH == 0) = 0;
                end
            end
            
            if idx ~= 1
                % recalculate voltage dependent load
                s = spu.*(aPQ + aI.*sqrt(E) + aZ.*E) - 1j*cappu + wpu;
                for k2 = nnode:-1:2
                    inline = nodes.inmat(k2);
                    outlines = nodes.outmat(nodes.outmat(:,k2) ~= 0,k2);
                    P(:,inline) = real(s(:,k2)) + sum(P(:,outlines),2);
                    Q(:,inline) = imag(s(:,k2)) + sum(Q(:,outlines),2);
                end
                P(lines.PH == 0) = 0; Q(lines.PH == 0) = 0;
            end
            % check if idx is a junction node
            if(isempty(find(jn == idx)))
                % then this is not a junction node, move up path
%                 disp 'move on';
                cn = idx;
%                 idx = find(FM(cn,:) == -1);
                idx = lines.TXnum(nodes.inmat(cn));
            else
                % this is a junction node and we stop this branch
%                 disp 'end path, junction reached';
                flag = false;
            end
%             if idx == 1
%                 flag = false;
%             end
%             pause
        end
        % completed branches from terminal nodes to junction nodes
        % junction nodes are now terminal branches
    end
    
    Etest = E(:,2);
    Dtest = D(:,2);
    iter = iter+1;
        
end

%     clc
%     disp('Junction')
%     E, D, P, Q
%     pause

V = sqrt(E).*exp(1j*pi/180*D);
V(nodes.PH == 0) = 0;

Stx = P + 1j*Q; Stx(lines.PH == 0) = 0;
Srx = Stx - Lmn; Srx(lines.PH == 0) = 0;

sV = spu.*(aPQ + aI.*sqrt(E) + aZ.*E) - 1j*cappu + wpu;

Vrot = [V(1,:); V(2,:).*exp(1j*120*pi/180); V(3,:).*exp(1j*240*pi/180)];

FBS.E = E;
FBS.D = D;
FBS.P = P;
FBS.Q = Q;

FBS.V = V;
FBS.Vrot = Vrot;
FBS.Stx = Stx;
FBS.Srx = Srx;
FBS.sV = sV;
FBS.iter = iter;

% FBS.I = I;
% FBS.Irot = Irot;
% FBS.Inode = Inode;

% for k1 = 1:1:nline
%     I(:,k1) = lines.FYpu(:,:,k1)*(V(:,lines.TXnum(k1)) - V(:,lines.RXnum(k1)));
% end
% I(lines.PH == 0) = 0;
% Irot = [I(1,:); I(2,:).*exp(1j*120*pi/180); I(3,:).*exp(1j*240*pi/180)];

% Currents consumed at each node
% Inode = conj(sV./V); Inode(nodes.PH==0) = 0;

% Power delivered to each node
% for k1=1:1:nline
%     Stx(:,k1) = V(:,lines.TXnum(k1)).*conj(I(:,k1));
%     Srx(:,k1) = V(:,lines.RXnum(k1)).*conj(I(:,k1));
% end
% Srx
% Stx

end