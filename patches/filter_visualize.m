function filter_visualize(filters)
for i = 1:numel(filters.half_patch_size)
  if(strcmp(filters.type{i}, 'local'))
    image = [];
    for j = 1:size(filters.filter{i},2)
      filt = filters.filter{i}(:,j);
      filt = filt/norm(filt);
      filter_shape = reshape(filt, 2*filters.half_patch_size(i)+1, 2*filters.half_patch_size(i)+1, []);
      im_row = [];
      for c = 1:size(filter_shape,3)
        im_row = [im_row, filter_shape(:,:,c)];
      end
      image = [image; im_row];
    end
    imwrite(image+.5, ['filter_image_' num2str(i) '.png']);
  end
end
end
