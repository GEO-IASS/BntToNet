%% Setup
addpath matpower4.1/
addpath matpower4.1/extras
addpath bnt/
addpath(genpathKPM('bnt'))

%% INPUT FROM MATPOWER

 a = case4;
% a = case4gs;
% a = case14;
% a = case118;
% a = case300;
bus = a.bus;
gen = a.gen;
branch = a.branch;

% if using DC
opt = mpoption('PF_DC', 1);
data = runpf(a, opt);

%% Covert to BNT
fprintf('Converting power flow to BNT format...\n')
tic
bnet = convertToBNT(data, struct('verbose',true));
bntToNet(bnet, 'out.net');
toc

%% Initialize alg
fprintf('Initializing junction tree...\n')
tic
engine = jtree_inf_engine(bnet);
toc
 
%% Do inference
fprintf('Running inference...\n')
tic
N = length(bnet.parents);
% w/ no clamping
evidence = cell(1,N);
[engine, ll] = enter_evidence(engine, evidence);
marg = marginal_nodes(engine, 1);
toc
