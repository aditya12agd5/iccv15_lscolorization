function [filter_vals] = apply_filter(filter_patches, filters, params)
warning('Depreciated');
filter_vals = [];
idx = 1;
for i = 1:params.filters.nums
  half_patch_width = params.filters.half_patch_width(i);
  window_size = (2*half_patch_width+1)^2;
  assert(window_size == numel(filters{i}));

  for j = 1:size(params.filters.locations{i},1)
    for c = 1:params.nchannels
      filter_vals(:,end+1) = filter_patches(:,idx:idx+window_size-1)*filters{i};
      idx = idx+window_size;
    end
  end

end

end
