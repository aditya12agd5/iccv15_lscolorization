function ColorizeMain(colorizetrain, colorizetest, output_dir, stage)
%% E.g. Usage: ColorizeMain('./data/beach_small/Train/', './data/beach_small/Test/', './data/ColorizeOut/');

addpath(genpath('../lscolorization/'));

%% sanity check on args
if nargin < 3
    disp('Error: Insufficient arguments');
    return;
elseif nargin == 3
    stage = 1;
else
    if(~exist([output_dir '/model.mat'], 'file') ||...
       ~exist([output_dir '/params.mat'], 'file') ||...
       ~exist([output_dir '/exemplar_hist.mat'], 'file') )
            disp('Error: No trained model to load');
            return;
    end  
    stage = 2;
end


mkdir(output_dir);

params = struct(...
'proj_dimensions', 12, ... % Number of dimensions to project down to.
'weps', .25, ...
'maxtries', 8, ... % Num trees to train.
'iter_per_tree', 2, ... % number of iteration steps to take per tree.
'features_per_tree', 10000, ...
'maxdepth', 60, ...
'maxsplitloc', 250, ... % Determines the number of split values to try for a single projection.
'minnode', 100, ... % How many pixels should there be in our leafs at a minimum.
'maxnodegradient', 200, ... % How many pixels should we use to compute gradients in a leaf (for speed, set to features_per_tree, to use all).
'maxgradientstep', 0.001, ...
'stratified_sample', false, ...
'fperi', 5000, ... % Features to be sampled per train image
'pyramid_debug', false, ...
'filter_half_patch', [1 3], ...
'num_filters', [6 6], ...
'filter_half_patch_bar_spot', [1 3 5],...
'filter_bar_orient', [0, pi/3, 2*pi/3], ...
'ms_bandwidth', 0.3 ...
);

if(stage == 1) %Train and Test
  [forest_model, params, exemplar_hist] = train_model(params, colorizetrain);
  params.initial_model = [];
  params.image_spec = [];
  save([output_dir '/model.mat'], 'forest_model');
  save([output_dir '/params.mat'], 'params');
  save([output_dir '/exemplar_hist.mat'], 'exemplar_hist');
  model_inference(params, forest_model, colorizetest, output_dir, 'test', exemplar_hist);
  %%uncomment below to run inference on train images
  %model_inference(params, forest_model, colorizetrain, output_dir, 'train', exemplar_hist);
elseif(stage == 2) %Only Test
  load([output_dir '/model.mat']);
  load([output_dir '/params.mat']);
  load([output_dir '/exemplar_hist.mat']);
  model_inference(params, forest_model, colorizetest, output_dir, 'test', exemplar_hist);
end

end

function [forest_model, params, exemplar_hist] = train_model(params, datatrain)
[~, ~, ~, ~, ~, ~, gegfmat, gegsmat, gegtmat, gegrmat, filters, params, exemplar_hist] = makefblockcolorize(datatrain, params);
forest_model = learnfromfblockgvsplit(gegfmat, gegsmat, gegtmat, filters, params);
end

function [err_out] = model_inference(params, forest, data_dir, output_dir, type, exemplar_hist)

d = dir([data_dir]);
d = d(arrayfun(@(x)x.name(1), d)~='.');

num_to_recon = numel(d);
for j = 1:num_to_recon
   [dcell, icell, ~] = readdatacolorize(data_dir, params, d(j).name);
   [pcell, pcell_split, ~, ~, ~, ~, reconstructions, reconstructions_raw, ~, params] = makeretfeatsgrawinscolorize(icell, dcell, params, false);

   image_spec_local = params.image_spec;
   for k = 1:numel(params.image_spec)
     image_spec_local(k).func = @(x) params.image_spec(k).func(x, dcell{1,1});
   end

   [image_recon, pred_reconstruction] = reconallpixforest(dcell, icell, pcell, pcell_split,...
     forest, params, 1, [output_dir '/' type '_' num2str(j) '_'], image_spec_local, exemplar_hist);
end

err_out = [];
end

