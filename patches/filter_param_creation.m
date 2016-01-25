function filter_params = filter_param_creation(filter_sizes, number_per, init_funcs)
filter_params.nums = numel(filter_sizes);

filter_params.half_patch_width = filter_sizes;

if(numel(number_per) == 1)
  number_per = ones(numel(filter_sizes),1)*number_per;
end

% Check that we don't have a mismatch in size now.
assert(numel(number_per) == numel(filter_sizes));

for i = 1:numel(filter_sizes)
  % Only make one filter centered at the pixel (for debugging)
  %filter_params.locations{i} = round(zeros(1, 2)*filter_sizes(i));
  filter_params.locations{i} = zeros(1, 2);%round(randn(number_per(i), 2)*filter_sizes(i));
  g = init_funcs{i}(1+2*filter_sizes(i));
  filter_params.initial{i} = g(:);
end
end
