function [filters, params] = make_reconstruction_features(rcell, masks, params)

if(isfield(params, 'filter_override') && ~isempty(params.filter_override))
  keyboard;
  filters.filter = params.filter_override;
  filters.type = params.filter_override_type;
  filters.half_patch_size = params.filter_override_half_patch_size;
  return;
end

filters.filter = {};
filters.type = {};
filters.half_patch_size = [];
for i = 1:numel(params.filter_half_patch)
  filter_sub_size = get_reconstruction_filters(rcell, masks, params.filter_half_patch(i), params.num_filters(i));

  for k = 1:params.nchannels
    for j = 1:params.num_filters(i)
      filters.filter = {filters.filter{:}, filter_sub_size(:, (k-1)*params.nchannels+j)};
      filters.type = {filters.type{:}, 'local'};
      filters.half_patch_size(end+1) = params.filter_half_patch(i);
    end
  end
end

filters.filter = {eye(params.nchannels), filters.filter{:}};
filters.type{end+1} = 'local';
filters.half_patch_size = [0, filters.half_patch_size];

for scale = params.filter_half_patch_bar_spot
  for orientation = params.filter_bar_orient
    bar = make_gabor_two(scale, 2, 2/scale, orientation);
    bar = kron(eye(params.nchannels), bar(:));
    filters.filter = {filters.filter{:}, bar};
    filters.type = {filters.type{:}, 'local'};
    filters.half_patch_size(end+1) = scale;
  end

  f_log = fspecial('log', 2*scale+1, (2*scale+1)/8);
  f_log = kron(eye(params.nchannels), f_log(:));
  filters.filter{end+1} = f_log;
  filters.type = {filters.type{:}, 'local'};
  filters.half_patch_size(end+1) = scale;

  gauss = fspecial('gaussian', 2*scale+1, (2*scale+1)/4);
  gauss = kron(eye(params.nchannels), gauss(:));
  filters.filter{end+1} = gauss;
  filters.type = {filters.type{:}, 'local'};
  filters.half_patch_size(end+1) = scale;
end

if isfield(params, 'filter_sim_half_patch')
  for scale = params.filter_sim_half_patch
    filters.filter{end+1} = ones(2*scale+1, 1)/(2*scale+1);
    filters.type{end+1} = 'sparse';
    filters.half_patch_size(end+1) = scale;
  end
end

%% uncomment to generate filter visualization images
%filter_visualize(filters);

end
