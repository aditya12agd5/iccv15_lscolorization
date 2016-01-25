function [tree] = learntree(forest_model, egfmat, egsmat, egtmat, params, egweights)

proj_dimensions = params.proj_dimensions;
reconstruction_size = params.reconstruction_size;
weps = params.weps;
nfeats = params.nfeats;


tfn = size(egfmat,1);
troot = node([],[], [1:tfn]', [], [], []);
tfnref=1:tfn;

if(nargin < 6)
  egweights = ones(size(egfmat,1), 1);
end

fscales = forest_model.fscales;

[egam, egbv] = get_M_and_B(forest_model, egsmat, params);

[~, gstars, egbs, ~] = evallrfun(egfmat, egtmat, egam, egbv, params, egweights);

% Compute gradients
wgbv = -gstars+(1/(1-weps))*(egbs-weps*gstars);
wgbm = zeros(tfn, nfeats*proj_dimensions);

% Last nzero entries get 0's for B and b.

%QUESTION: should this be egsmat?
for i = 1:proj_dimensions-params.nzeros
    wgbm(:, nfeats*(i-1)+1:nfeats*i) = ((wgbv(:, i))*ones(1, nfeats)).*(egfmat);
end

wgam=zeros(tfn, reconstruction_size*proj_dimensions);
for i=1:reconstruction_size
    wgam(:, proj_dimensions*(i-1)+1:proj_dimensions*i)=((1+weps)*gstars-egbs+...
        (1/(1-weps))*(egbs-weps*gstars)*(-weps)).*...
        ((egtmat(:, i))*ones(1, proj_dimensions));
    % i.e. zero gradient for negative values
    % this is the gradient with respect to the A at each location,
    % and is correct, by earlier codes
    % but we also need the jacobian of the A at the location wrt D
    % fortunately, this is zero for all but nceninds
end

eggv = [wgam, wgbm, wgbv];

% Make a tree
tree = makesimpletreeongv2(troot, egsmat, eggv, 0, fscales, params);

leafs = tree.flatten();

display('**** Tree ***');
% Initialize the tree we just built.
for leaf_idx = 1:numel(leafs)
    leafs(leaf_idx).amat = zeros(1, proj_dimensions*reconstruction_size);
    leafs(leaf_idx).bvmat = zeros(1, proj_dimensions*nfeats + proj_dimensions);
    leafs(leaf_idx).alpha = 1;
end

for iter = 1:params.iter_per_tree
  for leaf_idx = 1:numel(leafs)
    scost = evallrfun(egfmat, egtmat, egam, egbv, params, egweights);

    % now do an update on egmm, egbv
    weights = leafs(leaf_idx).weights;
    ceninds = tfnref(weights>0)';

    if numel(ceninds) > 3
      [egam(ceninds,:), egbv(ceninds,:), alpha, gam, gbv] =...
          updateambvnormzerosng(egam(ceninds,:), egbv(ceninds, :),...
              egfmat(ceninds, :), egtmat(ceninds, :), ...
              weights(ceninds, :), params, egweights(ceninds,:));

      leafs(leaf_idx).amat = leafs(leaf_idx).amat - gam * alpha;
      leafs(leaf_idx).bvmat = leafs(leaf_idx).bvmat - gbv * alpha;
      leafs(leaf_idx).alpha = 1;

      [wcost, gstars, egbs, vals] = evallrfun(egfmat, egtmat, egam, egbv, params, egweights);

      % Debug print
      display([num2str(wcost) ' ' num2str(scost)]);
      if wcost>scost+scost*1e-5 % went uphill by more than a trivial amount
        warning(['Step too large. ' num2str(wcost) ' from ' num2str(scost) ' alpha ' num2str(alpha)]);% this should never happen
      end
     % else
     %   disp(['Step good. ' num2str(wcost) ' from ' num2str(scost) ' alpha ' num2str(alpha)]);
     % end

    else
      leafs(leaf_idx).amat = zeros(1, proj_dimensions*reconstruction_size);
      leafs(leaf_idx).bvmat = zeros(1, proj_dimensions*nfeats + proj_dimensions);
      leafs(leaf_idx).alpha = 0;
      warning('Leaf too small.');
    end
  end
end

% Save memory, we don't need weights or inda after here.
for leaf_idx = 1:numel(leafs)
  leafs(leaf_idx).inda = [];
  leafs(leaf_idx).weights = [];
end
end
