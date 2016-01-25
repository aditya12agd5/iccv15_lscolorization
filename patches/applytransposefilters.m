function [im] = applytransposefilters(filtered_im, pcell, filters)
nchannels = size(filters.filter{1},1)/(2*filters.half_patch_size(1)+1)^2;
im = zeros(size(filtered_im,1), size(filtered_im,2), nchannels);

count = 1;
for j = 1:numel(filters.half_patch_size)
  if(strcmp('local', filters.type{j}))
    filter_width = 2*filters.half_patch_size(j)+1;
    filter_numel = filter_width^2;
    for k = 1:size(filters.filter{j},2)
      for c = 1:nchannels
        filter = reshape(filters.filter{j}(1+(c-1)*filter_numel:c*filter_numel,k), filter_width, filter_width);
        filter_flip = doubleflip(filter);
        im(:,:,c) = im(:,:,c) + imfilter(filtered_im(:,:,count), filter_flip);
      end
      count = count + 1;
    end
  elseif(strcmp('sparse', filters.type{j}))
    for c = 1:nchannels
      im(:,:,c) = im(:,:,c) + applytransposesimilarityfilter(pcell{1}, filtered_im(:,:,count), struct('half_patch_size', filters.half_patch_size(j), 'filter', filters.filter{j}));
      count = count + 1;
    end
  end
end

end
