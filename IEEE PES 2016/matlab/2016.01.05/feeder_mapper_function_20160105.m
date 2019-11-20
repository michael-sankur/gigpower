function [feeder, loads, controllers] = feeder_mapper_function_20160105(fp,fn1,fn2,fn3,fn4,fn5)

% Michael Sankur - msankur@berkeley.edu
% 2016.01.05

%% Load all feeder data

% Read from csv file, input first 4 columns as strings (general case)
% Change file path as needed

% Base value data
fid = fopen([fp fn1]);
DATA1 = textscan(fid,'%*s%f','delimiter',',','headerlines',1);
fclose(fid);

% Feeder topology
fid = fopen([fp fn2]);
DATA2 = textscan(fid,'%s%s%f%s','delimiter',',','headerlines',1);
fclose(fid);

% Feeder segment impedances
fid = fopen([fp fn3]);
DATA3 = textscan(fid,'%s%*s%f%f','delimiter',',','headerlines',1);
fclose(fid);

% Feeder loads
fid = fopen([fp fn4]);
DATA4 = textscan(fid,'%s%*s%s%f%f%f%f%f%f','delimiter',',','headerlines',1);
fclose(fid);

% Feeder controllers
fid = fopen([fp fn5]);
DATA5 = textscan(fid,'%s%s%f%f%f','delimiter',',','headerlines',1);
fclose(fid);

%% Parse Data

%% Assign base values

feeder.Vbase = DATA1{1}(1)*1e3;
feeder.Sbase = DATA1{1}(2)*1e3;
feeder.Ibase = feeder.Sbase/feeder.Vbase;
feeder.Zbase = feeder.Vbase/feeder.Ibase;

%% Find root node and tail nodes

nodeA = DATA2{1};
nodeB = DATA2{2};
seglength = DATA2{3};
segconf = DATA2{4};

for k1=1:length(nodeA)
    % loop through rows of A, find entry in A is not found in B. rn
    % will be the row of [A B] with the root node in A.
    if isempty(strmatch(nodeA(k1),nodeB,'exact')) == 1
        feeder.rn = k1;
        break;
    end
end
feeder.rn; % index of root node in A
feeder.RN = nodeA(feeder.rn); % root node names

feeder.tn = [];
for k1=1:length(nodeB)
    % loop through rows of B, find where entries in B are not found in A.
    % tn will be the rows of [A B] with tail nodes in B. this will only
    % work properly on if nodes only have one (1) parent
    if isempty(strmatch(nodeB(k1),nodeA,'exact')) == 1
        feeder.tn(end+1) = k1;
    end
end
feeder.tn; % indices of tail nodes in B
feeder.TN = nodeB(feeder.tn); % tail node names

%% Assign numbers to nodes

feeder.nodelist = [feeder.RN]; % list of nodes, with root node as 1
feeder.SL = [0]; % segment lengths
feeder.conflist = {'X'};
% loop through all relationships
for k1 = 1:length(nodeB)    
    if isempty(strmatch(nodeB(k1),feeder.nodelist,'exact')) == 1
        feeder.nodelist(end+1) = nodeB(k1);
        feeder.SL(end+1) = seglength(k1); % add segment length to SL
        feeder.conflist(end+1) = segconf(k1);
    end
end

n = length(feeder.nodelist);
feeder.n = n;

%% Create feeder relationship matrix and matrix of paths

feeder.FM = zeros(n,n); % matrix of node relationships
% iterate down rows of [A B] (node relationships) and add relationships to
% FM, and segment lengths to S
for k1 = 1:length(nodeA)
    % upstream (nA) and downstream (nB) node for relationship
    nA = nodeA(k1);
    nB = nodeB(k1);
    % find where in nodelist nA and nB are
    idxA = strmatch(nodeA(k1),feeder.nodelist,'exact');
    idxB = strmatch(nodeB(k1),feeder.nodelist,'exact');
    
    % add relationship to FM
    feeder.FM(idxA,idxB) = 1;
    feeder.FM(idxB,idxA) = -1;
        
end

% lists of node numbers
feeder.tn = []; % list of tail nodes
feeder.pn = []; % list of pass through nodes
feeder.jn = []; % list of junction nodes
for k1 = 1:n
    if sum(feeder.FM(k1,:)) == -1
        feeder.tn(end+1) = k1;      
    elseif sum(feeder.FM(k1,:)) == 0
        feeder.pn(end+1) = k1;
    elseif sum(feeder.FM(k1,:)) >= 1
        feeder.jn(end+1) = k1;
    end    
end

% matrix of all paths from feeder head to tail
% paths that are shorter than the max have zeros added onto them
% use find(k,:) function to obtain kth path
feeder.paths = [];
for k1 = 2:n
    if sum(feeder.FM(k1,:)) == -1
        temppath = k1;
        while temppath(end) ~= 1
            temppath(end+1) = find(feeder.FM(temppath(end),:) == -1);
        end
        temppath = fliplr(temppath);
        if length(temppath) > size(feeder.paths,2)
            feeder.paths = [feeder.paths zeros(size(feeder.paths,1), length(temppath)-size(feeder.paths,2)); temppath];
        elseif length(temppath) == size(feeder.paths,2)
            feeder.paths = [feeder.paths; temppath];
        elseif length(temppath) < size(feeder.paths,2)
            feeder.paths = [feeder.paths; temppath zeros(1,size(feeder.paths,2)-length(temppath))];
        end
    end    
end

%% Create phase and impedance matrices

confname = DATA3{1};
segR = DATA3{2};
segX = DATA3{3};

% matrices of dimension 3x3xn with feeder impedance
feeder.FZpl = zeros(3,3,n); % feeder segment impedance per unit length
feeder.FZ = zeros(3,3,n); % feeder segment impedance
feeder.FZpu = zeros(3,3,n); % feeder segment impdenace per unit
for k1 = 2:n
    idx = strmatch(feeder.conflist(k1),confname,'exact');
    idx = idx(1);
    feeder.FZpl(:,:,k1) = [segR(idx) + j*segX(idx), segR(idx+1) + j*segX(idx+1), segR(idx+2) + j*segX(idx+2);
        segR(idx+1) + j*segX(idx+1), segR(idx+3) + j*segX(idx+3), segR(idx+4) + j*segX(idx+4);
        segR(idx+2) + j*segX(idx+2), segR(idx+4) + j*segX(idx+4), segR(idx+5) + j*segX(idx+5)];
    feeder.FZ(:,:,k1) = feeder.FZpl(:,:,k1)*feeder.SL(k1)/5280;
end
feeder.FZpu = feeder.FZ/feeder.Zbase;

% Create matrix of phases that exist at nodes, rows indicate phase,
% columns indicate node
% 0 = no phase exists at node, 1 = phase exists at node
feeder.PH = zeros(3,n);
feeder.PH(:,1) = [1 1 1]';
for k1 = 2:n
    for ph = 1:3
        if feeder.FZ(ph,ph,k1) ~= 0
            feeder.PH(ph,k1) = 1;
        end
    end
end

clear confname segR segX

%% Create Load and Capacitor Matrices

loads.s = zeros(3,n);
loads.cap = zeros(3,n);
for k1 = 1:length(DATA4{1});
    k2 = strmatch(DATA4{1}(k1),feeder.nodelist,'exact');
    if strmatch(DATA4{2}(k1),'load','exact')
        loads.s(:,k2) = 1e3*[DATA4{3}(k1) + j*DATA4{4}(k1);
            DATA4{5}(k1) + j*DATA4{6}(k1);
            DATA4{7}(k1) + j*DATA4{8}(k1)];
    elseif strmatch(DATA4{2}(k1),'cap','exact')
        loads.cap(:,k2) = 1e3*[j*DATA4{4}(k1);
            j*DATA4{6}(k1);
            j*DATA4{8}(k1)];
    end
end

loads.spu = loads.s/feeder.Sbase;
loads.cappu = loads.cap/feeder.Sbase;

%% Create controllers and parameters

abc = {'a','b','c'};

controllers.cn = zeros(3,n); % matrix of existing DER at phase/node
controllers.wmax = zeros(3,n); % matrix of maximum DER power
controllers.rr = zeros(3,n); % matrix of DER rate rate percent of max magnitude
for k1 = 1:length(DATA5{1});
    k2 = strmatch(DATA5{1}(k1),feeder.nodelist,'exact');
    ph = strmatch(DATA5{2}(k1),abc,'exact');
    controllers.cn(ph,k2) = 1;
    controllers.cnstate(ph,k2) = DATA5{3}(k1);
    controllers.wmax(ph,k2) = DATA5{4}(k1);
    controllers.rr(ph,k2) = DATA5{5}(k1);
end

controllers.cnodes = find(sum(controllers.cn)~=0);

controllers.rr = controllers.rr/100.*controllers.wmax;

controllers.wmaxpu = controllers.wmax/feeder.Sbase;
controllers.rrpu = controllers.rr/feeder.Sbase;

%% Cleanup

clear fp fn1 fn2 fn3 fn4 fn5
clear ph k1 k2 nA nB temppath  A B C D E F G H I DATA1 DATA2 DATA3 DATA4 DATA5

end