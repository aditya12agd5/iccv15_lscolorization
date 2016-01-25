function [filters_out] = pad_filters(filters, filter_half_width)

filter_count = 1;
for i = 1:numel(filters.half_patch_size)
  for j = 1:size(filters.filter{i}, 2)
    filter_image = reshape(filters.filter{i}(:,j), 2*filters.half_patch_size(i)+1, 2*filters.half_patch_size(i)+1, []);
    filter_image_padded = padarray(filter_image, [filter_half_width-filters.half_patch_size(i), filter_half_width-filters.half_patch_size(i)], 'replicate', 'both');
    filters_out.filter{1}(:,filter_count) = filter_image_padded(:);
    filter_count = filter_count + 1;
  end
end
filters_out.half_patch_size(1) = filter_half_width;
end
