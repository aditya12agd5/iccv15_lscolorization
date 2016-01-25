function [images, recon_raw, recon_indep] = reconallpixforest(dcell, icell, pcell, pcell_split, forest_model, params, ntorecon, image_prefix, image_spec, exemplar_hist)
% images is the image spec applied to the reconstruction
% recon_raw is the reconstruction reshaped as an image without any function applied.

% Get the params local variables
nfeats = params.nfeats;
proj_dimensions = params.proj_dimensions;
reconstruction_size = params.reconstruction_size;

images = cell(ntorecon, length(image_spec));

for imptr=1:ntorecon
    fprintf(1, 'Reconstructing ****************** %d\n', imptr);
    [nums, numt, ~] = size(icell{imptr,1});

    [b_img2, b_img2_masked, b_img_transpose_filt, nwmat, mask, wmbv] = recon_matrix_construction(icell{imptr,1}, dcell{imptr,end}, pcell{imptr,1}, pcell_split{imptr, 1}, forest_model, params);

    % Solve for min objective.
    res = pcg(@(u)applymat_filter(nwmat, pcell{imptr,:}, forest_model.filters, nums, numt, 1e-5, u), b_img2, 1e-5, 30);

    [~, true_hist, ~] = obtainbmyrmghist(icell{imptr, 2});
    
    %last histogram is the true histogram for the given image 
    exemplar_hist{end+1} = int32(round(true_hist));
 
	im_hist_res = cell(1, numel(exemplar_hist));
	im_hist_res_ext = cell(1, numel(exemplar_hist));

	for nHist = [1:numel(exemplar_hist)]
		[res_after_hist, corresp, bDistInit] = histogramOptimize(res, nwmat, b_img2, wmbv, pcell{imptr, :}, forest_model.filters, nums, numt, 1e-5, exemplar_hist{nHist}, params.nchannels, params.ms_bandwidth, image_prefix);
   	 	im_hist_res{nHist} = reshape(reshape(res_after_hist, params.nchannels, [])', nums, numt, []);
		if(nHist == numel(exemplar_hist))
			im_hist_res_ext{nHist} = ['-gt_hist_' num2str(nHist) '-']; 
		else
			im_hist_res_ext{nHist} = ['-hist_' num2str(nHist) '-'];
		end 
	end

    im_res = reshape(reshape(res, params.nchannels, [])', nums, numt, []);
    recon_raw{imptr} = im_res;
  
    im_b_indep = [];
    for nchannel = 1:params.nchannels
    	im_b_indep(:, :, nchannel) = reshape(wmbv(nchannel:params.proj_dimensions:end), nums, numt); 
    end 

    for i = 1:length(image_spec)
      indep_img = image_spec(i).func(b_img_transpose_filt(:,:,image_spec(i).sample));
      recon_indep{imptr, i} = indep_img.*repmat(dcell{imptr,end}, [1,1, size(indep_img,3)]);
      images{imptr, i} = image_spec(i).func(recon_raw{imptr}(:,:,image_spec(i).sample)).*repmat(dcell{imptr,end}, [1,1,size(indep_img,3)]);

      im_hist_res_rgb = cell(1, numel(exemplar_hist));
      
      for nHist = [1:numel(exemplar_hist)]
	      im_hist_res_rgb{nHist} = image_spec(i).func(im_hist_res{nHist}).*repmat(dcell{imptr,end}, [1,1,size(indep_img,3)]);
      end

      im_b_indep_rgb = image_spec(i).func(im_b_indep).*repmat(dcell{imptr,end}, [1,1,size(indep_img,3)]);

      if(isfield(image_spec(i), 'name') && ~strcmp(image_spec(i).name, ''))
        imwrite(images{imptr, i}, [image_prefix image_spec(i).name '-' num2str(imptr) '.tif']);
        imwrite(recon_indep{imptr, i}, [image_prefix image_spec(i).name '-indep-' num2str(imptr) '.tif']);

	for nHist = [1:numel(exemplar_hist)]
	        imwrite(im_hist_res_rgb{nHist}, [image_prefix image_spec(i).name im_hist_res_ext{nHist} num2str(imptr) '.tif']);
	end

        imwrite(im_b_indep_rgb, [image_prefix image_spec(i).name '-b_indep-' num2str(imptr) '.tif']);
      end

    end
end

end

