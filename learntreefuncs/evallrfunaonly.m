function [ecost, gstars, vals] = evallrfunaonly(egtmat, egam, egbs, params, weights)

weps = params.weps;

% Params
proj_dimensions = params.proj_dimensions;
reconstruction_size = params.reconstruction_size;

%gstars = egtmat;
gstars=zeros(size(egtmat, 1), proj_dimensions);
for i=1:proj_dimensions
    gstars(:, i)=sum(egtmat.*egam(:, 1+(i-1)*reconstruction_size:i*reconstruction_size), 2);
end
vals = sum(gstars.^2, 2)*(1+weps)*(1/2)-sum(gstars.*egbs, 2)+(1/(1-weps))*(1/2)*(sum((egbs-weps*gstars).^2, 2));
assert(all(vals >= 0));

if(nargin < 5)
  warning('unweighted evallrfunaonly');
  ecost=sum(vals);
else
  ecost = sum(vals.*weights);
end
