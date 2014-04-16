%% erikreed@cmu.edu
clear all
DO_INFERENCE = false;
DC = 0;

if ~DC
  fprintf('Running AC power\n');
else
  fprintf('Running DC power\n');
end

cases = dir('matpower4.1/case*');
runtimeConvert = zeros(1,length(cases));
runtimeInit = zeros(1,length(cases));
runtimeInf = zeros(1,length(cases));
stats = cell(1,length(cases));

timestr = '%.2f seconds.\n';

outdir = 'out';
fprintf('Saving results to ./%s', outdir)
mkdir(outdir)

for i=1:length(cases)
  name = cases(i).name;
  if ~isempty(strfind(name, 'caseformat'))
    fprintf('Skipping %s\n', name)
    continue
  end
  fprintf('Running: %s\n', name)
  
  run(name);
  a = ans;

  opt = mpoption('PF_DC', DC);
  
  [~, data] = evalc('runpf(a, opt)');
  
  try
    % Covert to BNT
    fprintf('Converting power flow to BNT format... ')
    tic
    bnet = convertToBNT(data);
    stats{i} = bnet.stats;
    
    runtimeConvert(i) = toc;
    fprintf(timestr, runtimeConvert(i))
    memory
    
    basename = strrep(name, '.m', '');
    dotname = sprintf('%s/%s.dot', outdir, basename);
    fprintf('Saving %s\n', dotname);
    % TODO: color-coding script
    graph_to_dot(bnet.dag, 'leftright', 1, 'filename', ...
        dotname','node_label', bnet.names)
    
    netname = sprintf('%s/%s.net', outdir, basename);
    bntToNet(bnet, netname);
    
    if DO_INFERENCE
      % Initialize alg
      fprintf('Initializing junction tree... ')
      tic
      engine = jtree_inf_engine(bnet);
      runtimeInit(i) = toc;
      fprintf(timestr, runtimeInit(i))
      memory
      
      % Do inference
      fprintf('Running inference... ')
      tic
      N = length(bnet.parents);
      % w/ no clamping
      evidence = cell(1,N);
      [engine, ll] = enter_evidence(engine, evidence);
      marg = marginal_nodes(engine, 1);
      runtimeInf(i) = toc;
      fprintf(timestr, runtimeInf(i))
      memory
      
      fprintf('Completed %s\n', name)
    end
    
  catch ME
    fclose('all');
    idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match');
    if isempty(idSegLast)
       rethrow(ME)
    end
    switch idSegLast{1}
      case 'nomem'
        fprintf('Out of memory error on %s!\n', name)
      otherwise
        rethrow(ME)
    end
  end
  clear bnet data marg evidence N a engine name
end
fclose('all');
fprintf('Done!\n\n')

%% print stats in table format
if DO_INFERENCE
  fprintf('Name, Buses, Branches, Loads, Generators, Nodes, Convert Time, Compile Time, Inference Time\n')
  for i=1:length(cases)
    if ~isempty(stats{i})
      fprintf('%s, %d, %d, %d, %d, %d, %.2f, %.2f, %.2f\n', cases(i).name, ...
        stats{i}.buses, stats{i}.branches, stats{i}.loads, ...
        stats{i}.generators, stats{i}.nodes, runtimeConvert(i), runtimeInit(i), ...
        runtimeInf(i));
    end
  end
  
  fprintf('\n\nFailed bnets:\n');
  for i=1:length(cases)
    if isempty(stats{i})
      fprintf('%s\n', cases(i).name);
    end
  end
end
