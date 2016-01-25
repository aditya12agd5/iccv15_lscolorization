function [dcell, icell, exemplar_hist] = readdatacolorize(colorize_root, params, image_dir_name)

count = 1;

dnames=dir(colorize_root);
if(nargin == 3)
	dnames = struct('name', image_dir_name);
end
dnames = dnames(arrayfun(@(x)x.name(1), dnames)~='.');

count_hist = 0;

for i = 1:numel(dnames)

  eps = 0.001; %resetting epsilon

  full_images_dir = [colorize_root '/' dnames(i).name]

  color_image = imresize(im2double(imread([full_images_dir '/color.png'])), [256 NaN]);
  gray_image = sum(color_image, 3)./3;
  bmy_rmg_img =  converttobmyrmg(color_image); 
  
  mask = im2double(ones(size(gray_image)));
  mask = double(mask > 0.1);

  dcell{count, 1} = gray_image;
  dcell{count ,2} = bmy_rmg_img;
  dcell{count, 3} = mask;
  
  icell{count, 1} = gray_image;
  icell{count, 2} = bmy_rmg_img;
  icell{count, 3} = log(gray_image + eps);
  icell{count, 4} = log(bmy_rmg_img + eps);

  [~, curr_hist, ~] = obtainbmyrmghist(bmy_rmg_img);

  if(i == 1)
	mean_hist = curr_hist;
	count_hist = count_hist + 1;
  else
  	mean_hist = mean_hist + curr_hist;
	count_hist = count_hist + 1;
  end
  
  avg_bmy_rmg_image = double(zeros(size(bmy_rmg_img)));

  store_match_img = {};

  numMatch = 9;
  for j = [1:numMatch] 
  	match_color_image = im2double(imread([full_images_dir '/match' int2str(j) '.jpg']));
  	match_color_image = imresize(match_color_image, [size(gray_image, 1) size(gray_image,2)]);
  	
	%(@aditya): to deal with noisy match images in all-scenes dataset
  	if(size(match_color_image, 3) == 1)
  		match_color_image = repmat(match_color_image, [1 1 3]);	
  	end
  	
  	match_bmy_rmg_img = converttobmyrmg(match_color_image);
  	store_match_img{j} = match_bmy_rmg_img;
  	avg_bmy_rmg_image = avg_bmy_rmg_image + match_bmy_rmg_img; 
  end
  avg_bmy_rmg_image = avg_bmy_rmg_image./numMatch;

  icell{count, 5} = (avg_bmy_rmg_image);
  
  var = double(zeros(size(bmy_rmg_img)));
  for j = [1:numMatch]
  	var = var+((store_match_img{j} - avg_bmy_rmg_image).^2);
  end
  var = var./numMatch;
  icell{count, 6} = sqrt(var);
  
  count = count + 1;
end

%%uses only mean histogram, can add other histograms here
exemplar_hist = {};
exemplar_hist{end+1} = int32(round(mean_hist./count_hist)); 

end
