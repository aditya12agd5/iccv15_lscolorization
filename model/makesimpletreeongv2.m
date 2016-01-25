function nnode = makesimpletreeongv2(current_node, egsmat, eggv, depth, fscales, params)

% TODO(jrock): Pretty confident that the weights really aren't necessary in
% training especially since we are making hard decissions about which leaf
% nodes go into here, and then just reweighting after the hard split is decided.

spflag=0;
if (depth<params.maxdepth)&&(current_node.nindices>params.minnode)
  [left, right, spflag]=spliterrorongmean(current_node, egsmat, eggv, fscales, params);
end

if spflag==1
  ltree = makesimpletreeongv2(left, egsmat, eggv, depth+1, fscales, params);
  rtree = makesimpletreeongv2(right, egsmat, eggv, depth+1, fscales, params);
  nnode = node(ltree, rtree, [], ltree.cumsplitprojection(:, 1:end-1), ltree.cumsplitval(1:end-1), ltree.cumsplitsign(1:end-1));
else
  nnode=current_node;
  % it is a leaf
  nnode.leaf=1;
  % make the weights corresponding to this leaf
  weights = nnode.getWeights(egsmat, fscales);
  nnode.weights = weights;
end
