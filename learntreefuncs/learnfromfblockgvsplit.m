function [forest_model] = learnfromfblockgvsplit(egfmat, egsmat, egtmat, filters, params)
% Params
proj_dimensions = params.proj_dimensions;
nfeats = params.nfeats;
nfeats_split = params.nfeats_split;

% This is the number of features to try when creating splits,
% and add this to params.

if(~isfield('ntos', params))
%  params.ntos = min(100, floor(sqrt(nfeats)) + 1);
   params.ntos = min(100, floor(sqrt(nfeats_split)) + 1);
   if(nfeats_split < 50)
 	  params.ntos = nfeats_split;
   end
end

[tree0, params] = params.initial_model(params);
forest = {tree0};

% fscales rescales the features for splitting so that the sigmoid is roughly the same
% size for any feature.

fscales_size = 25;
fscales=zeros(nfeats_split, 1);
for i=1:nfeats_split
  if(max(egsmat(:,i)) == min(egsmat(:,i)))
    fscales(i) = fscales_size;
  else
    fscales(i)=fscales_size./(max(egsmat(:, i))-min(egsmat(:, i))); % so it's 100 steps from smallest to largest
  end
end


for i = 1:params.maxtries

  forest_model = struct('forest', {forest}, 'fscales', fscales, 'filters', filters);

  %%IGNORE stratified_sample case for now
  if params.stratified_sample
    sample = randsample(size(egfmat,1), min(params.features_per_tree*25, size(egfmat,1)));
    segfmat = egfmat(sample, :);
    segtmat = egtmat(sample, :);

    display('create stratified sample');
    [egam, egbv] = get_M_and_B(forest_model, segfmat, params);

    [~, ~, ~, numerators] = evallrfun(segfmat, segtmat, egam, egbv, params);

    [~, idx] = sort(numerators);
    sample = [];
    prev_split = 1;
    split_ind = log([2,3,4,5,6])/log(6);

    for split = floor(split_ind*numel(idx))
      disp(['cost at split: ' num2str(numerators(idx(split)))]);
      sample = [sample; randsample(idx(prev_split:split-1), min(1/(numel(split_ind)) * params.features_per_tree, split-prev_split))];
      prev_split = split;
    end

    numel(sample)

    egweights = 1./numerators(sample);
    egweights = numel(egweights)*egweights/norm(egweights);
    new_tree = learntree(forest_model, segfmat(sample, :), segtmat(sample, :), params, egweights);
  else
    display('create random sample');
    sample = randsample(size(egfmat,1), min(params.features_per_tree, size(egfmat,1)));
   
    params_new = params;
    params_new.maxdepth = floor(i/params.maxtries)*params.maxdepth;
    params_new.maxgradientstep = max(1, ceil(5-5*(i/100)));
    params_new.iter_per_tree = max(1, ceil(5-5*(i/100)));
    new_tree = learntree(forest_model, egfmat(sample, :), egsmat(sample, :), egtmat(sample, :), params_new);

  end
  forest{end + 1} = new_tree;
end

forest_model = struct('forest', {forest}, 'fscales', fscales, 'filters', filters);

end
