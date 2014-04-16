function bnet = convertToBNT(data, opts)

LOW_MEMORY = true;

% -- defaults --
verbose = false;
discretize = false;
num_discrete_states = 1;
% -------------

if isfield(opts, 'verbose')
    verbose = opts.verbose;
end
if isfield(opts, 'discretize')
    discretize = opts.discretize;
end
if isfield(opts, 'num_discrete_states')
    num_discrete_states = opts.num_discrete_states;
end

assert(mod(num_discrete_states, 2) == 1)

%Convert data from matpower to BNT

%% Fix bad IDs
USE_SEQUENTIAL_IDS = false;

ids = sort(unique(vertcat(data.bus(:,1), data.bus(:,2), data.gen(:,1), ...
    data.branch(:,1), data.branch(:,2))));
numIDs = length(ids);
idMap = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
idMapInv = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
for i=1:numIDs
    idMap(ids(i)) = i;
    idMapInv(i) = ids(i);
end

if USE_SEQUENTIAL_IDS
    for i=1:numIDs
        idMapInv(i) = i;
    end
end

for i=1:length(data.bus(:,1))
    data.bus(i,1) = idMap(data.bus(i,1));
end

for i=1:length(data.branch(:,1))
    data.branch(i,1) = idMap(data.branch(i,1));
    data.branch(i,2) = idMap(data.branch(i,2));
end

for i=1:length(data.gen(:,1))
    data.gen(i,1) = idMap(data.gen(i,1));
end

bus = data.bus;
gen = data.gen;
branch = data.branch;

%% INPUT

B = size(data.bus,1);
G = size(data.gen,1);
BR = size(data.branch,1);
L = sum(data.bus(:, 3)~=0);
loadIds = find(data.bus(:, 3) ~= 0);
N = B+G+BR*3+L;

if verbose
    fprintf('Number of branches: %d\n', BR)
    fprintf('Number of generators: %d\n', G)
    fprintf('Number of buses: %d\n', B)
    fprintf('Number of loads: %d\n', L)
    fprintf('Number of BN nodes: %d\n', N)
end

if LOW_MEMORY && N > 3000
    fprintf('Low memory and many nodes! Terminating...\n');
    throw(MException('ResultChk:nomem', 'Low memory and many nodes!'));
end

%% Initialize DAG

%node 1~G:gen
%node G+1~G+BR:branch
%node G+GR+1~N:bus
dag = zeros(N,N);

%Gen->Bus
for i = 1:G
    dag(i,G+BR*3+L+gen(i,1)) = 1;
end

%Load->Bus
loadsAdded = 0;
for l = 1:B
    if data.bus(l, 3) ~= 0
        loadsAdded = loadsAdded + 1;
        dag(loadsAdded+G,G+BR*3+L+l) = 1;
    end
end
assert(loadsAdded == L);

%LineLoss->LineNet
for ll = 1:BR
    dag(ll+G+L+BR,G+L+BR*2+ll) = 1;
end

%BR->Bus or LineNet
for k = 1:BR
    if data.branch(k,14)>0
        dag(k+L+G,G+BR*3+L+branch(k,1)) = 1;
        dag(k+L+G,G+L+BR*2+k) = 1;
    else
        dag(k+L+G,G+L+BR*2+k) = 1;
        dag(k+L+G,G+BR*3+L+branch(k,2))=1;
    end
end

%LineNet->bus
for k = 1:BR
    if data.branch(k,14)>0
        dag(G+L+BR*2+k,G+BR*3+L+branch(k,2)) = 1;
    else
        dag(G+L+BR*2+k,G+BR*3+L+branch(k,1)) = 1;
    end
end

%% Add CPTs
names = {};
ns = ones(1,N) * num_discrete_states;

discrete_nodes = [];
if discretize
    discrete_nodes = 1:N;
end

bnet = mk_bnet(dag, ns, 'discrete', discrete_nodes);
if discretize
    display 'TODO: Discretized BN CPDs are all random at the moment'
    return
end
for x = 1:G
    if ~discretize
        bnet.CPD{x} = gaussian_CPD(bnet, x,'mean', data.gen(x,2));
    else
        error('unimplemented')
    end
    names{length(names)+1} = sprintf('Gen_%d', x);
end

for x = 1:L
    bnet.CPD{x+G} = gaussian_CPD(bnet,x+G,'mean',data.bus(loadIds(x),3));
    % To get ordering corresponding to load number rather than load ID,
    % use data.bus(loadIds(x))
    names{length(names)+1} = sprintf('Load_%d', x);
end

for x = 1:BR
    bnet.CPD{x+G+L} = gaussian_CPD(bnet,x+G+L,'mean', ...
        max(abs(data.branch(x,[14 16]))));
    names{length(names)+1} = sprintf('Line_%d', x);
end

for x = 1:BR
    % Note: this will always be 0 when running DC
    bnet.CPD{x+G+L+BR} = gaussian_CPD(bnet,x+G+L+BR,'mean', ...
        abs(data.branch(x,14)+data.branch(x,16)));
    names{length(names)+1} = sprintf('Line_%d_Loss', x);
end

for x = 1:BR
    parent = cell2mat(bnet.parents(x+G+L+2*BR));
    mean_net(x) = abs(data.branch(parent(1)-G-L,14));
    assert(length(parent) == 2);
    bnet.CPD{x+G+L+2*BR} = gaussian_CPD(bnet,x+G+L+BR*2,'mean',mean_net(x), ...
        'weights', [1, -1]);
    names{length(names)+1} = sprintf('Line_%d_Net', x);
end

% +gen, -line, -load, +linenet
for x = 1:B
    p = cell2mat(bnet.parents(x+G+L+3*BR));
    n = size(p,2);
    node_mean = 0;
    weights = [];
    for m = 1:n
        if p(m)<G+1
            node_mean = node_mean + data.gen(p(m),2);
            weights(m) = 1;
        elseif p(m)<G+L+1
            node_mean = node_mean - data.bus(p(m)-G,3);
            weights(m) = -1;
        elseif p(m)<G+L+BR+1
            node_mean = node_mean - max(abs(data.branch(p(m)-G-L,[14 16])));
            weights(m) = -1;
        elseif p(m)>G+L+2*BR
            node_mean = node_mean + mean_net(p(m)-G-L-2*BR);
            weights(m) = 1;
        end
    end
    
    bnet.CPD{x+G+L+3*BR} = gaussian_CPD(bnet,x+G+L+BR*3,'mean',node_mean, ...
        'weights', weights);
    names{length(names)+1} = sprintf('Bus_%d', x);
end

bnet.names = names;

bnet.stats = struct();
bnet.stats.branches = BR;
bnet.stats.buses = B;
bnet.stats.loads = L;
bnet.stats.generators = G;
bnet.stats.nodes = N;

end
