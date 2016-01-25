function [b_img2, b_img2_masked, b_img_transpose_filt, nwmat, mask, wmbv] = recon_matrix_construction(image, mask, pyr, pyr_split, forest_model, params)

proj_dimensions = params.proj_dimensions;
reconstruction_size = params.reconstruction_size;
nfeats = params.nfeats;

nums = size(image,1);
numt = size(image,2);

% Get the 2D pointer as well as the linear pointer.
[iptrs, jptrs] = find(mask);
ptrs = find(mask);


gegfmat = pyr_to_feature_vec(pyr);
gegsmat = pyr_to_feature_vec(pyr_split);

[gegam, gegbv] = get_M_and_B(forest_model, gegsmat(ptrs,:), params);

% I now have the am and the bv, bvmat at each pixel
tstep=1; % always
nic=size(1:tstep:nums, 2);
njc=size(1:tstep:numt, 2);
%chicken!
ntiles=nums*numt;
nnz=size(ptrs, 1);

wmbv=zeros(proj_dimensions*ntiles, 1);
mmi=zeros(proj_dimensions*reconstruction_size*nnz, 1);
mmj=zeros(proj_dimensions*reconstruction_size*nnz, 1);
wmvals=zeros(proj_dimensions*reconstruction_size*nnz, 1);
smwp=1;

for i=1:nnz
  vptr=ptrs(i);
  viptr=iptrs(i);
  vjptr=jptrs(i);


  tnumber=(vjptr-1)*(nums)+(viptr);
  op=proj_dimensions*(tnumber-1)+1;
  wb1=(reshape(gegbv(i, 1:proj_dimensions*nfeats), nfeats, proj_dimensions))'*gegfmat(vptr, :)';
  wb2=reshape(gegbv(i, proj_dimensions*nfeats+1:end), proj_dimensions, 1);

  wm=(reshape(gegam(i, :), reconstruction_size, proj_dimensions))';

  assert(~any(isnan(wb2+wb1)) && ~any(isnan(wm(:))));

  wmbv(op:op+proj_dimensions-1)=wb2+wb1;
  for j=1:proj_dimensions
      mmi(smwp:smwp+reconstruction_size-1)=(op+j-1)*ones(reconstruction_size, 1);
      mmj(smwp:smwp+reconstruction_size-1)=((tnumber-1)*reconstruction_size+1)+[0:reconstruction_size-1]';
      wmvals(smwp:smwp+reconstruction_size-1)=wm(j, :);
      smwp=smwp+reconstruction_size;
  end
end

%Sparse matrix construction
nwmat = sparse(mmi, mmj, wmvals, proj_dimensions*ntiles, reconstruction_size*ntiles);
%[wcols, patches, weave] = wovenwcolsmaker(nums, numt, params);

b_img = reshape(reshape((nwmat'*wmbv)', reconstruction_size, [])', nums, numt, reconstruction_size);
b_img_transpose_filt = applytransposefilters(b_img, pyr, forest_model.filters);
b_img2 = reshape(reshape(b_img_transpose_filt, [], params.nchannels)', [], 1);
b_img2_masked = b_img_transpose_filt(repmat(mask, [1,1,params.nchannels]));

end

