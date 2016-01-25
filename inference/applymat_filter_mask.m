function vec=applymat(amat, filters, nums, numt, deps, vec, mask, default_val)

im = default_val;
im(mask) = vec;
%im = reshape(reshape(vec, [], nums*numt)', nums, numt, []);
%im = reshape(vec, nums, numt, []);

filtered_im = applyfilters(im, filters);
im_vec = reshape(reshape(filtered_im,[],size(filtered_im,3))',[],1);
v1 = amat*im_vec;
v2 = amat'*v1;
im_reshape = reshape(reshape(v2,size(filtered_im,3),[])', nums, numt, size(filtered_im,3));
im_transpose_filters = applytransposefilters(im_reshape, filters);

vec = im_transpose_filters(mask) + deps*vec;
%vec = reshape(reshape(im_transpose_filters,[],size(im_transpose_filters,3))',[],1) + deps * vec;

%vec = amat'* (amat * vec) + deps * vec;
%vec = wcols'*(nwmat'*(nwmat*(wcols*vec)))+deps*vec;
% TODO(jrock): not sure what to make of this.
% (w'm'mw)*v + d*v
