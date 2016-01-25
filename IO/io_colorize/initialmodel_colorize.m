function [nnode, params] = initialmodel_colorize(params)
iam=zeros(params.proj_dimensions, params.reconstruction_size);

%since we have two channel image
setrows = 2;

% This selects the current pixel.  4 channels match to Reflectance RGB and Shading.
iam(1:setrows, 1:setrows) = eye(setrows);
iam(setrows+1:params.proj_dimensions,:) = .01*randn(params.proj_dimensions-setrows, params.reconstruction_size);

ibv=zeros(params.proj_dimensions, 1); % no constant term
ibvf=zeros(params.proj_dimensions, params.nfeats);


ibvf(1, 10)=1; % average bmy estimate
ibvf(2, 11)=1; % average rmg estimate

nnode = node([], [], [], zeros(params.nfeats_split, 1), [-1], [1]);
nnode.leaf = true;

% A matrix reshaped as a vector.
startam=reshape(iam', 1, []);

% B matrix reshaped and a vector and then append b.
startbv=[reshape(ibvf', 1, []), reshape(ibv, 1, [])];

nnode.amat = startam;
nnode.bvmat = startbv;
nnode.alpha = 1;

% Params
params.nzeros = 0;
params.setrows = setrows;

