function [FBS] = FBS_TX(feeder,nodes,lines,configs,loads,caps,controllers,Vref)

% Michael Sankur - msankur@berkeley.edu
% 2016.06.03
% Adapted from function originally by Dan Arnold

% Output is the FBS struct.
% FBS.V is the 3xn matrix of complex pu voltages phasors
% FBS.I is the 3xn matrix of complex pu current phasors
% FBS.Inode is the 3xn matrix of complex pu current phasors delivered to
% the node
% FBS.S is the 3xn matrix of complex pu power phasors delivered to the node
% from an upstream node
% FBS.sV is the 3xn matrix of complex pu power phasors consumed by the node
% FBS.iter is the number of iterations

% Feeder parameters
% n = feeder.n;
% FM = feeder.FM;
% PH = feeder.PH;
% FZpu = feeder.FZpu;
% jn = feeder.jn;
% tn = feeder.tn;

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

% Setup voltages and currents
% V = [Vref Vref zeros(3,nnode-2)];
V = Vref*ones(1,nnode);
I = zeros(3,nline);
Inode = zeros(3,nnode);
S = zeros(3,nnode);

% Initialize iteration
iter = 0;
tol = abs(Vref(2))*1e-9;
Vtest = zeros(3,1);

while(max(abs(Vtest-Vref)) >= tol)
    
    V(:,2) = Vref;
    % sweep forward, calculate voltages
    for k1 = 2:1:nnode
        % find children of current node
        lineidx = nodes.outmat(nodes.outmat(:,k1) ~= 0,k1);
        for k2 = 1:1:length(lineidx)
            kline = lineidx(k2);
            knode = lines.RXnum(kline);
%             kline = nodes.outmat(k2,k1);
            % calculate voltage at child node
            V(:,knode) = V(:,k1) - FZpu(:,:,kline)*I(:,kline);
            % zero out voltage for nonexistent phases at all nodes
            V(nodes.PH == 0) = 0;
        end
        clear idx;
    end
            
    % recalculate voltage dependent load
    s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V).^2)) - 1j*cappu + wpu;
        
    % TERMINAL NODES
%     disp 'Terminal Nodes'
    % loop through terminal nodes
    for k1 = length(tn):-1:1
%         disp('NEW TERMINAL NODE')
        cn = tn(k1);
        kline = nodes.inmat(cn);
        % recalculate voltage dependent load
        s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V).^2)) - 1j*cappu + wpu;
        % calculate segment current
%         s(:,cn)
%         V(:,cn)
%         s(:,cn)./V(:,cn)
%         conj(s(:,cn)./V(:,cn))
%         Iline = conj(s(:,cn)./V(:,cn))
        I(:,kline) = conj(s(:,cn)./V(:,cn));
        % zero out current for nonexistent phases at all nodes
        I(lines.PH==0) = 0;
                
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
                   V(ph,idx) = V(ph,cn) + FZpu(ph,:,kline)*I(:,kline);
                end
            end
            % zero out voltage for nonexistent phases at all nodes
            V(nodes.PH==0) = 0;
            
            % recalculate voltage dependent load
            s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V).^2)) - 1j*cappu + wpu;
            % calculate current in parent segment
            I(:,nodes.inmat(idx)) = I(:,kline) + conj(s(:,idx)./V(:,idx));
            % zero out current for nonexistent phases at all nodes
            I(lines.PH==0) = 0;
            
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
    s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V).^2)) - 1j*cappu + wpu;
        
    % JUNCTION NODES
    % start from bottom-most junction node and repeat until a junction
    % node is hit, stop if junction node is 1
    % loop through junction nodes
%     disp 'Junction Nodes'
    for k1 = length(jn):-1:2
%         disp('NEW JUNCTION NODE')
        cn = jn(k1);
        kline = nodes.inmat(cn);
        % calculate segment current
        % for all nodes except 1, more than 1 segment is connected
        % calculate node current
        Itot(:,1) = conj(s(:,cn)./V(:,cn));
        % find children of current node
%         idx = find(FM(cn,:) == 1);
%         idx = nodes.outmat(nodes.outmat(:,cn) ~= 0,cn);
        lineidx = nodes.outmat(nodes.outmat(:,cn) ~= 0,cn);
        for k2 = 1:1:length(lineidx)
            % add current of childrens' segment (if any) to current segment
            % current
            Itot = Itot + I(:,lineidx(k2));
        end
        I(:,kline) = Itot;
        % zero out current for nonexistent phases at all nodes
        I(lines.PH==0) = 0;
        
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
                    V(ph,idx) = V(ph,cn) + FZpu(ph,:,kline)*I(:,kline);
                end
            end
            % zero out voltage for nonexistent phases at all nodes
            V(nodes.PH==0) = 0;
            
            if idx ~= 1
                % recalculate voltage dependent load
                s = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V).^2)) - 1j*cappu + wpu;
                % calculate current in parent segment
                I(:,nodes.inmat(idx)) = I(:,kline) + conj(s(:,idx)./V(:,idx));
                % zero out current for nonexistent phases at all nodes
                I(lines.PH==0) = 0;
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
                
    Vtest = V(:,2);
    iter = iter+1;
        
end

% Power delivered to each node
for k1=1:1:nline
    Stx(:,k1) = V(:,lines.TXnum(k1)).*conj(I(:,k1));
    Srx(:,k1) = V(:,lines.RXnum(k1)).*conj(I(:,k1));
end

% Voltage dependent complex loads
sV = spu.*(aPQ + aI.*(abs(V)) + aZ.*(abs(V).^2)) - 1j*cappu + wpu;

% Currents consumed at each node
Inode = conj(sV./V); Inode(nodes.PH==0) = 0;

Vrot = [V(1,:); V(2,:).*exp(1j*120*pi/180); V(3,:).*exp(1j*240*pi/180)];
Irot = [I(1,:); I(2,:).*exp(1j*120*pi/180); I(3,:).*exp(1j*240*pi/180)];

FBS.V = V;
FBS.Vrot = Vrot;
FBS.I = I;
FBS.Irot = Irot;
FBS.Inode = Inode;
FBS.Stx = Stx;
FBS.Srx = Srx;
FBS.sV = sV;
FBS.iter = iter;

end