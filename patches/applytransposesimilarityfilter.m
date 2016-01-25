function [im] = applytransposesimilarityfilter(im, filtered_im, filter)
[similarity, M] = getsimilarity(im, size(filtered_im,3), 2*filter.half_patch_size + 1, filter.filter);

im = M'*filtered_im(:);
im = reshape(im, size(filtered_im));
end

