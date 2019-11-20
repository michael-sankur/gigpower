function [FBS] = FBS_3phase_fun_20160105(feeder,loads,control,V0)

% Michael Sankur - msankur@berkeley.edu
% 2016.01.05
% Adapted from function originally by Dan Arnold

% feeder parameters
n = feeder.n;
FM = feeder.FM;
PH = feeder.PH;
FZpu = feeder.FZpu;
jn = feeder.jn;
tn = feeder.tn;

% load parameters
spu = loads.spu;
A0 = loads.A0;
A1 = loads.A1;
cappu = loads.cappu;

% controller parameters
wpu = control.wpu;

% setup voltages and currents
V = [V0 zeros(3,n-1)];
I = zeros(3,n);
Inode = zeros(3,n);
S = zeros(3,n);

% Initialize iteration
iter = 0;
tol = abs(V0(1))*0.0001;
Vtest = zeros(3,1);

while(max(abs(Vtest-V0)) >= tol)
    
    V(:,1) = V0;
    % sweep forward, calculate voltages
    for k1 = 1:1:n
        % find children of current node
        idx = find(FM(k1,:)==1);        
        for k2 = 1:1:length(idx)
            % calculate voltage at child node
            V(:,idx(k2)) = V(:,k1) - FZpu(:,:,idx(k2))*I(:,idx(k2));
            % zero out voltaget for nonexistent phases at all nodes
            V(PH == 0) = 0;
        end
        clear idx;
    end
    
    % recalculate voltage dependent load
    s = spu.*(A0 + A1.*(abs(V).^2)) - cappu + wpu;

    % TERMINAL NODES
%     disp 'Terminal Nodes'
    % loop through terminal nodes
    for k1 = length(tn):-1:1
        cn = tn(k1);        
        % recalculate voltage dependent load
        s = spu.*(A0 + A1.*(abs(V).^2)) - cappu + wpu;
        % calculate segment current
        I(:,cn) = conj(s(:,cn)./V(:,cn));
        % zero out current for nonexistent phases at all nodes
        I(PH==0) = 0;
        
        % move up the feeder until a junction node is hit        
        % assume each node only has 1 parent
        flag = true;
        idx = find(FM(cn,:) == -1);
        while(flag)
            % calculate voltage at parent node only for phases existing at
            % current node. doing so for all phases will set a phase to
            % an erroneous value if it exists at the parent and not at the
            % current node
            for ph = 1:3
                if PH(ph,cn) == 1
                   V(ph,idx) = V(ph,cn) + FZpu(ph,:,cn)*I(:,cn);
                end
            end
            % zero out voltaget for nonexistent phases at all nodes
            V(PH==0) = 0;
            
            % recalculate voltage dependent load
            s = spu.*(A0 + A1.*(abs(V).^2)) - cappu + wpu;
            % calculate current in parent segment
            I(:,idx) = I(:,cn) + conj(s(:,idx)./V(:,idx));
            % zero out current for nonexistent phases at all nodes
            I(PH==0) = 0;
            
            % check if idx is a junction node
            if(isempty(find(jn == idx)))
                % then this is not a junction node, move up path
%                 disp 'move on';
                cn = idx;
                idx = find(FM(cn,:) == -1);
            else
                % this is a junction node and we stop this branch
%                 disp 'end path, junction reached';
                flag = false;
            end            
        end
        % completed branches from terminal nodes to junction nodes
        % junction nodes are now terminal branches
    end
    
    % recalculate voltage dependent load
    s = spu.*(A0 + A1.*(abs(V).^2)) - cappu + wpu;
    
    % JUNCTION NODES
    % start from bottom-most junction node and repeat until a junction
    % node is hit, stop if junction node is 1
    % loop through junction nodes
%     disp 'Junction Nodes'
    for k1 = length(jn):-1:2
        cn = jn(k1);
        % calculate segment current
        % for all nodes except 1, more than 1 segment is connected
        % calculate node current
        Itot(:,1) = conj(s(:,cn)./V(:,cn));
        % find children of current node
        idx = find(FM(cn,:) == 1);
        for k2 = 1:1:length(idx)
            % add current of childrens' segment (if any) to current segment
            % current
            Itot = Itot + I(:,idx(k2));
        end
        I(:,cn) = Itot;
        % zero out current for nonexistent phases at all nodes
        I(PH==0) = 0;
        
        % move up the feeder until a junction node is hit
        flag = true;
        idx = find(FM(cn,:) == -1);        
        while(flag)
            % calculate voltage at parent node only for phases existing at
            % current node. doing so for all phases will set a phase to
            % an erroneous value if it exists at the parent and not at the
            % current node
            for ph = 1:3
                if PH(ph,cn) == 1
                   V(ph,idx) = V(ph,cn) + FZpu(ph,:,cn)*I(:,cn);
                end
            end
            % zero out voltaget for nonexistent phases at all nodes
            V(PH==0) = 0;
            
            % recalculate voltage dependent load
            s = spu.*(A0 + A1.*(abs(V).^2)) - cappu + wpu;
            % calculate current in parent segment
            I(:,idx) = I(:,cn) + conj(s(:,idx)./V(:,idx));
            % zero out current for nonexistent phases at all nodes
            I(PH==0) = 0;
            
            % check if idx is a junction node
            if(isempty(find(jn == idx)))
                % then this is not a junction node, move up path
%                 disp 'move to'
                cn = idx;
                idx = find(FM(cn,:) == -1);
            else
                % this is a junction node and we stop this branch
%                 disp 'end path, junction reached'
                flag = false;
            end
        end
        % completed branches from terminal nodes to junction nodes
        % junction nodes are now terminal branches
    end
        
    Vtest = V(:,1);
    iter = iter+1;
    
%     pause
    
end

% power delivered to each node
for k1=1:1:n
    S(:,k1) = V(:,k1).*conj(I(:,k1));
end

% voltage dependent loads
sV = spu.*(A0 + A1.*(abs(V).^2)) - cappu + wpu;

% currents consumed at each node
Inode = sV./V; Inode(PH==0) = 0;


FBS.V = V;
FBS.I = I;
FBS.Inode = Inode;
FBS.S = S;
FBS.sV = sV;
FBS.iter = iter;

end