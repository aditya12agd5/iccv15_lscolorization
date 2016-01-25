function col_im = im2col_multichannel(image, params)
% Pads so that size(col_im, 2) == size(image,1)*size(image,2)

image_pad = padarray(image, [params.half_patch_width, params.half_patch_width]);

col_im = [];
for i = 1:size(image_pad,3)
  col_im = [col_im; im2col(image_pad(:,:,i), [params.winsize, params.winsize])];
end

end
