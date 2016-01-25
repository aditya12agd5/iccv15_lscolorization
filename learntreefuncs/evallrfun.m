function [ecost, gstars, egbs, vals] = evallrfun(egfmat, egtmat, egam, egbv, params, weights)

% From Params.
proj_dimensions = params.proj_dimensions;
reconstruction_size = params.reconstruction_size;
weps = params.weps;

nfeats = params.nfeats;

%gstars = egtmat;
gstars=zeros(size(egtmat, 1), proj_dimensions);
for i=1:proj_dimensions
    gstars(:, i)=sum(egtmat.*egam(:, 1+(i-1)*reconstruction_size:i*reconstruction_size), 2);
end

egbs=egbv(:, nfeats*proj_dimensions+1:end);
for i=1:proj_dimensions
    egbs(:, i)=egbs(:, i)+sum(egfmat.*egbv(:, 1+(i-1)*nfeats:i*nfeats), 2);
end

vals=sum(gstars.^2, 2)*(1+weps)*(1/2)-sum(gstars.*egbs, 2)+(1/(1-weps))*(1/2)*(sum((egbs-weps*gstars).^2, 2));
assert(all(vals >= 0));

if(nargin < 6)
  warning('unweighted evallrfun');
  ecost = sum(vals);
else
  ecost = sum(weights.*vals);
end

assert(ecost >= 0);
