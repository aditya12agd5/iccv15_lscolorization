function [pcell, pcell_split, params] = makefeaturepyramidcolorize(icell, dcell, params, filters)

num_images = size(icell, 1);

% make pyramids for each image
pcell=cell(num_images, 1);
pcell_split = cell(num_images, 1);

for wim=1:num_images

    gray_image = icell{wim, 1}; %take the orig image, no log

    nums=size(gray_image, 1); numt=size(gray_image, 2); dx=1; dy=1;
    [axmat, aymat, axxmat, axymat, ayymat]=setupmatrices2fast(nums, numt, dx, dy);

    gray_image_vec = gray_image(:);
    lx=reshape(axmat*gray_image_vec, nums, numt, 1);
    ly=reshape(aymat*gray_image_vec, nums, numt, 1);

    % #total feature maps
    pyrcell=cell(1, 229);

    pyrcell(1:9)={gray_image, lx, ly, (lx.^2+ly.^2).^(1/2), ...
        imresize(imresize(gray_image, 1/2), [nums, numt]), ...
        imresize(imresize(gray_image, 1/4), [nums, numt]), ...
        imresize(imresize(gray_image, 1/8), [nums, numt]), ...
        imresize(imresize(gray_image, 1/16), [nums, numt]), ...
        imresize(imresize(gray_image, 1/32), [nums, numt])};

    bmy_rmg_img = icell{wim,5};
    pyrcell(10:11) = {bmy_rmg_img, icell{wim,6}}; %mean top k image and variance image
    
    %use normalized location as a feature
    [xCord yCord] = meshgrid([1:numt], [1:nums]);
    cords = [];
    cords(:, :, 1) = xCord./numt;
    cords(:, :, 2) = yCord./nums;
    pyrcell(12) = {cords};
   
    %use a global descriptor of gray image as feature
    global_feat_vec = reshape(imresize(gray_image, [5, 5]), 1, [])';
    for global_feat_dim = 1:numel(global_feat_vec)
    	pyrcell(12+global_feat_dim) = {repmat(global_feat_vec(global_feat_dim), nums, numt)};
    end

    %for lm filter responses
    lm_filters = make_lmfilters;
    startLm = 12+numel(global_feat_vec);
    for wth=1:size(lm_filters,3)
        conv_img = conv2(gray_image, lm_filters(:,:,wth), 'same');
        lm_response = (conv_img-min(conv_img(:))) ./ (max(conv_img(:)-min(conv_img(:))));

        lm_response_dx=reshape(axmat*lm_response(:), size(lm_response));
        lm_response_dy=reshape(aymat*lm_response(:), size(lm_response));

        %NO BLUR OF LM RESPONSE
        pyrcell(startLm+4*(wth-1)+1:startLm+4*wth) = {...
            lm_response, lm_response_dx, lm_response_dy,...
            (lm_response_dx.^2+lm_response_dy.^2).^(1/2)};
    end

    %use texture+global_image features as regression features
    pcell(wim, 1)={pyrcell(13:end)}; 
    
    %%use match color and grey intensity features as split features
    pcell_split(wim, 1) = {pyrcell(1:12)}; 
end

end
