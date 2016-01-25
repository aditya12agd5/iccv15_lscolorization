function [egam, egbv] = get_M_and_B(forest_model, egsmat, params)

proj_dimensions = params.proj_dimensions;
reconstruction_size = params.reconstruction_size;
nfeats = params.nfeats;
nfeats_split = params.nfeats_split;

forest = forest_model.forest;
fscales = forest_model.fscales;

npix = size(egsmat, 1);
%this is equal to size(egfmat, 1)

egam = zeros(npix, proj_dimensions * reconstruction_size);
egbv = zeros(npix, proj_dimensions * nfeats + proj_dimensions);

wind = 1;
for i = 1:numel(forest)
  tree = forest{i};
  leafs = tree.flatten();
  for j = 1:numel(leafs)
    if(numel(leafs(j).amat) == 0)
      continue;
    end

    %wv = leafs(j).getWeights(egfmat, fscales);
    wv = leafs(j).getWeights(egsmat, fscales);
    %wv_(egsmat,1)_1;
  
    sel = find(wv ~= 0);

    bvmat = leafs(j).bvmat*leafs(j).alpha;
    amat = leafs(j).amat*leafs(j).alpha;

    %tic;
   % egam_ = egam(sel,:);
   % egbv_ = egbv(sel,:);

   % parfor k = 1:numel(sel)
   %   egam_(k,:) = egam(k,:) + wv(sel(k))*amat;
   %   egbv_(k,:) = egbv_(k,:) + wv(sel(k))*bvmat;
   %   %egam(sel(k),:) = egam(sel(k),:) + wv(sel(k))*amat;
   %   %egbv(sel(k),:) = egbv(sel(k),:) + wv(sel(k))*bvmat;
   % end
   % egam(sel,:) = egam_;
   % egbv(sel,:) = egbv_;
    egam = egam + wv*amat;
    egbv = egbv + wv*bvmat;
   % egam(sel,:) = egam(sel,:) + wv(sel)*amat;
   % egbv(sel,:) = egbv(sel,:) + wv(sel)*bvmat;
    %egbv_add = tic;
   % for k = find(sel)
   %   egam(k,:) = egam(k,:) + wv(k)*amat;
   %   egbv(k,:) = egbv(k,:) + wv(k)*bvmat;
   % end
   % toc

    egam(sel,:) = orthonormalize(egam(sel,:), proj_dimensions, reconstruction_size);
    %egbv(sel,:)=egbv(sel,:)+(wv(sel)*leafs(j).bvmat)*leafs(j).alpha;
    %time_dist = toc(egbv_add);
    %display(['adding time took ' num2str(time_dist)]);
    % TODO(jrock): Remove for speed.
    %assert(~any(isnan(egbv(:))) && ~any(isnan(egam(:))));

    if (mod(wind, 100) == 0)
        fprintf(1, '%d\n', wind);
    end
    wind = wind + 1;
  end
end

end
