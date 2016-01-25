function [filters] = get_reconstruction_filters(reconstructions, masks, filter_half_width, num_filters)
filter_width = 2*filter_half_width + 1;


valid_patches = [];
for i = 1:length(reconstructions)
  col_im = im2col_multichannel(reconstructions{i}, struct('half_patch_width', filter_half_width, 'winsize', filter_width));
  mask_col = im2col_multichannel(masks{i},struct('half_patch_width', filter_half_width, 'winsize', filter_width));
  select = sum(mask_col) >= filter_width^2;
  sub_col_im = col_im(:,select);
  select_downsize = randsample(size(sub_col_im,2), min(size(sub_col_im,2),10000));
  valid_patches = [valid_patches, sub_col_im(:, select_downsize)];
end

% Select filters for a single channel.
filter_channel = {};
mau_channel = {};
pca_sample = randsample(size(valid_patches,2), 10000);
for i = 1:size(reconstructions{1},3)
  selection = valid_patches((i-1)*filter_width^2+1:i*filter_width^2,pca_sample)';
  [filters_pca, ~, variance, ~, ~, mu] = pca(selection);
  variance = sqrt(variance);
  full_variance = sum(variance);
  cum_variance = [variance(1)];
  for j = 2:numel(variance)
    cum_variance(j) = cum_variance(j-1) + variance(j);
  end
  cum_variance = cum_variance/full_variance;

  selected_variance = [1:num_filters];%(cum_variance < .95 & variance' > .01);
  filter_channel{i} = filters_pca(:,selected_variance)./repmat(variance(selected_variance)',size(filters_pca,1),1);
  mu_channel{i} = mu;
  values_channel{i} = selection * filter_channel{i};
end

filters = zeros(size(reconstructions{1},3)*filter_width^2, size(reconstructions{1},3)*size(filter_channel{i},2));

g = fspecial('gaussian', [filter_width, filter_width], filter_width/3);
count = 1;
for i = 1:size(reconstructions{1},3)
  for j = 1:size(filter_channel{i},2)
    filter_single = g .* reshape(filter_channel{i}(:,j), filter_width, filter_width);
    filters((i-1)*filter_width^2+1:i*filter_width^2,count) = filter_single(:);
    count = count + 1;
  end
end
return;

values_all_channel = cat(2, values_channel{:});

pca_all_channel = pca(values_all_channel);

filters_pca = zeros(size(reconstructions{1},3)*filter_width^2, num_filters);
for i = 1:num_filters
  count = 1;
  filter_combined = [];
  for j = 1:numel(filter_channel)
    filter = zeros(size(filter_channel{j},1),1);
    for k = 1:size(filter_channel{j},2)
      filter = filter + filter_channel{j}(:,k) * pca_all_channel(count,i);
      count = count + 1;
    end
    filter_combined = [filter_combined; filter];
  end
  filters_pca(:,i) = filter_combined;
end

filters = filters_pca;
g = fspecial('gaussian', [filter_width, filter_width], filter_width/3);
for i = 1:size(filters,2)
  filter_channel = reshape(filters(:,i), filter_width, filter_width, []);
  for c = 1:size(filter_channel,3)
    filter_channel(:,:,c) = g .* filter_channel(:,:,c);
  end
  filters(:,i) = filter_channel(:);
end
end
