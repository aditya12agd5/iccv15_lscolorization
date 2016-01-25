function [tegtmat, tegfmat, tegrmat, filters, params] = maketrainingreconstructions(rcell, pcell, masks, params)
[filters, params] = make_reconstruction_features(rcell, masks, params);

[tegtmat, tegfmat, tegrmat, params] = maketrainingdatafromfilters(rcell, pcell, masks, filters, params);
end
