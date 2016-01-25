function [similarity, M] = getsimilarity(sym_im, nchannels, number, filter)
X = reshape(sym_im, [], size(sym_im,3));
similarity = knnsearch(X, X, 'K', number)';
% [~, similarity] = pdist2(X, X, 'euclidean', 'Smallest', number);

if(nargout >= 2)
  I = repmat(1:size(similarity,2), number, 1);
  F = repmat(filter(:), size(similarity,2), 1);
  M = kron(eye(nchannels), sparse(I(:), similarity(:), F));
end

similarity = similarity';
end
