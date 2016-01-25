function [filtered_im] = applyfilters(im, pcell, filters)
filtered_im = [];


for j = 1:numel(filters.half_patch_size)
  if(strcmp('local', filters.type{j}))
    filter_width = 2*filters.half_patch_size(j)+1;
    filter_numel = filter_width^2;
    im_filtered = zeros(size(im));
    for k = 1:size(filters.filter{j},2)
      for c = 1:size(im,3)
        filter = reshape(filters.filter{j}(1+(c-1)*filter_numel:c*filter_numel,k), filter_width, filter_width);
        im_filtered(:,:,c) = imfilter(im(:,:,c), filter);
      end
      filtered_im = cat(3, filtered_im, sum(im_filtered,3));
    end
  elseif(strcmp('sparse', filters.type{j}))
    im_filtered = applysimilarityfilter(im, pcell{1}, struct('half_patch_size', filters.half_patch_size(j), 'filter', filters.filter{j}));
    filtered_im = cat(3, filtered_im, im_filtered);
  end
end

end
