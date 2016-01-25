function [filtered_im] = applysimilarityfilter(im, sym_im, filter)
[similarity, M] = getsimilarity(sym_im, size(im,3), 2*filter.half_patch_size + 1, filter.filter);

filtered_im = M*im(:);
filtered_im = reshape(filtered_im, size(im));
end
