%% erikreed@cmu.edu
function bntToNet(bnet, outname)

numNodes = length(bnet.parents);

assert(numNodes == length(bnet.CPD))
assert(numNodes == length(bnet.node_sizes))
assert(numNodes == length(bnet.parents))

if ~isempty(bnet.cnodes)
  fprintf('Continuous nodes detected. Output will be unreadable by SamIAm.\n');
end

if isempty(bnet.names)
  names = {};
  names.keys = arrayfun(@(x)(sprintf('N%d',x)), 1:numNodes, 'UniformOutput', false);
else
  if iscell(bnet.names)
    names = bnet.names;
  else
    names = struct(bnet.names);
    assert(length(names) == numNodes)
    assert(length(names.keys) == length(names.vals))
    assert(length(names.keys) == numNodes)
    names = names.keys;
  end
end

fprintf('Writing to %s\n', outname);
fid = fopen(outname, 'w');

% print nodes
for i=1:numNodes
  numStates = bnet.node_sizes(i); 
  nodeName = names{i};
  isContinuous = any(bnet.cnodes == i);
  
  if isContinuous
    nodePrefix = 'continuous node';
  else
    nodePrefix = 'node';
    assert(any(bnet.dnodes == i))
  end
  
  fprintf(fid, '%s %s\n{\n', nodePrefix, nodeName);
  fprintf(fid, '  label = "%s";\n', nodeName);
  if ~isContinuous
    fprintf(fid, '  states = (');
    for j=1:numStates
      fprintf(fid, '"s%d" ', j);
    end
    fprintf(fid, ');\n');
  end
  fprintf(fid, '}\n\n');
end

% print CPTs
for i=1:numNodes
  parents = bnet.parents{i};
  nodeName = names{i};

  if ~isempty(parents)
    fprintf(fid, 'potential (%s | ', nodeName);
    for j=1:length(parents)
      parentName = names{parents(j)};
      fprintf(fid, '%s ', parentName);
    end
  else
    fprintf(fid, 'potential (%s', nodeName);
  end
  fprintf(fid, ')\n{\n');
  fprintf(fid, '  data = ');
  

  cpt = struct(bnet.CPD{i});
  isContinuous = any(bnet.cnodes == i);
  
  if isContinuous
    fprintf(fid, 'normal ( ');
    if isempty(parents)
      fprintf(fid, '%f, %f', cpt.mean, abs(cpt.mean)/10);
    else
      for j=1:length(parents)
        parentName = names{parents(j)};
        assert(abs(cpt.weights(j)) == 1);
        prefix = '+';
        if cpt.weights(j) < 0
          prefix = '-';
        end
        fprintf(fid, '%s %s ', prefix, parentName);
      end
      fprintf(fid, ', 1');
    end
    fprintf(fid, ' )');
  else
    
    CPT = cpt.CPT;
    n = ndims(CPT);
    parents_size = size(CPT);
    parents_size = parents_size(1:end-1);
    
    last = ind2subv(parents_size, 1);
    currentParens = length(parents_size);
    
    for j=1:currentParens+1
      fprintf(fid, '(');
    end
    
    for j=1:prod(parents_size)
      parent_inst = ind2subv(parents_size, j);
      numDiff = length(parents_size);
      for k=1:length(parent_inst)
        if last(k) ~= parent_inst(k)
          break
        end
        numDiff = numDiff - 1;
      end
      for k=1:numDiff
        fprintf(fid, ')');
      end
      fprintf(fid, '\n');
      for k=1:numDiff
        fprintf(fid, '(');
      end

      index = num2cell([parent_inst 1]);
      index{n} = ':';
      fprintf(fid, '%s', num2str(CPT(index{:}), '%f '));
      
      last = parent_inst;
    end
    
    for j=1:length(parents_size)+1;
      fprintf(fid, ')');
    end
  end
  fprintf(fid, ';\n');
  fprintf(fid, '}\n\n');
end

fclose(fid);
