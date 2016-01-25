function [vec, v1]=applymat(amat, pcell, filters, nums, numt, deps, vec)

im = reshape(reshape(vec, [], nums*numt)', nums, numt, []);
%im = reshape(vec, nums, numt, []);

filtered_im = applyfilters(im, pcell, filters);
im_vec = reshape(reshape(filtered_im,[],size(filtered_im,3))',[],1);
v1 = amat*im_vec;
v2 = amat'*v1;
im_reshape = reshape(reshape(v2,size(filtered_im,3),[])', nums, numt, size(filtered_im,3));
im_transpose_filters = applytransposefilters(im_reshape, pcell, filters);

vec = reshape(reshape(im_transpose_filters,[],size(im_transpose_filters,3))',[],1) + deps * vec;

