function [tegtmat, tegfmat, tegsmat, tegrmat,  params] = maketrainingdatafromfilters(rcell, pcell, pcell_split, masks, filters, params)

filtered_cell = cell(size(rcell));
tegtmat = cell(size(rcell));
tegfmat = cell(size(rcell));
tegsmat = cell(size(rcell));
tegrmat = cell(size(rcell));

for i = 1:numel(rcell)
  filtered_cell{i} = applyfilters(rcell{i}, pcell{i}, filters);

  column_cell = reshape(filtered_cell{i},[],size(filtered_cell{i},3));
  mask_column_cell = im2col_multichannel(masks{i}, struct('half_patch_width', max(filters.half_patch_size), 'winsize', 2*max(filters.half_patch_size)+1));
  if(max(params.filter_half_patch) == 0)
    window_valid = mask_column_cell >= 1;
  else
    window_valid = sum(mask_column_cell) >= 1;
  end

  tegtmat{i} = column_cell(window_valid,:);
  tegfmat{i} = pyr_to_feature_vec(pcell{i}, window_valid);
  tegsmat{i} = pyr_to_feature_vec(pcell_split{i}, window_valid);

  assert(size(tegsmat{i}, 1) == size(tegfmat{i}, 1), ...
	'Each feature map must be same length.');
 
  tegrmat{i} = zeros(sum(window_valid),3);
  [wx, wy] = find(reshape(window_valid, size(masks{i})));
  tegrmat{i}(:, 1) = wx;
  tegrmat{i}(:, 2) = wy;
  tegrmat{i}(:, 3) = i;
end

params.reconstruction_size = size(tegtmat{1},2);
params.nfeats = size(tegfmat{1}, 2);
params.nfeats_split = size(tegsmat{1}, 2);

end

