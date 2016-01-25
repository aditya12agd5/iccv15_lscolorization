function [pcell, pcell_split, gegfmat, gegsmat, gegtmat, gegrmat, reconstructions, reconstructions_raw, filters, params] = makeretfeatsgrawinscolorize(icell, dcell, params, training)


if(nargin < 4)
  training = true;
end

num_indices = [2];
params.image_spec = createImageSpec(num_indices, @toRGB, {'recon_color'});
params.nchannels = sum(num_indices);

reconstructions = {};
reconstructions_raw = {};
masks = {};
%%

for imptr=1:size(icell,1) % from the first images
  
  image = icell{imptr, 2};
  mask = dcell{imptr, 3};

  reconstructions{imptr} = image;
  reconstructions_raw{imptr} = icell{imptr, 2}; %color image
  
  masks{imptr} = mask;
end

if training

  gegtmat = [];
  gegfmat = [];
  gegsmat = [];
  gegrmat = [];

  [filters, params] = make_reconstruction_features(reconstructions, masks, params);
  for imptr=1:size(icell,1)
    disp(['image ' num2str(imptr) ' of ' num2str(size(icell,1))]);
 
    [pcell, pcell_split, params] = makefeaturepyramidcolorize(icell(imptr,:), dcell(imptr,:), params);
  
    [tegtmat, tegfmat, tegsmat, tegrmat, params] = maketrainingdatafromfilters({reconstructions{imptr}}, pcell, pcell_split, {masks{imptr}}, filters, params);
   
    if(params.fperi>size(tegtmat{1},1))
      gegfmat = [gegfmat; tegfmat{1}];
      gegsmat = [gegsmat; tegsmat{1}];
      gegtmat = [gegtmat; tegtmat{1}];
      gegrmat = [gegrmat; tegrmat{1}];
    else
      selection = randsample(size(tegtmat{1},1), params.fperi);
      gegfmat = [gegfmat; tegfmat{1}(selection, :)];
      gegsmat = [gegsmat; tegsmat{1}(selection, :)];
      gegtmat = [gegtmat; tegtmat{1}(selection, :)];
      gegrmat = [gegrmat; tegrmat{1}(selection, :)];
    end
 end

else
  [pcell, pcell_split, params] = makefeaturepyramidcolorize(icell, dcell, params);

  gegtmat = [];
  gegfmat = [];
  gegsmat = [];
  gegrmat = [];
  filters = [];
end
end
