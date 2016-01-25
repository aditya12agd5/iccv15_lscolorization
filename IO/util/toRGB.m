function [im_rgb] = toRGB(im_bmy_rmg, im_gray)

im_rgb = double(zeros(size(im_gray,1), size(im_gray,2), 3));

%red color
im_rgb(:,:,1) = im_gray - ((0.33).*(im_gray.*im_bmy_rmg(:,:,1))) + (0.5.*(im_gray.*im_bmy_rmg(:,:,2)));

%green color
im_rgb(:,:,2) = im_gray - ((0.33).*(im_gray.*im_bmy_rmg(:,:,1))) - (0.5.*(im_gray.*im_bmy_rmg(:,:,2)));

%blue color
im_rgb(:,:,3) = im_gray + ((0.66).*(im_gray.*im_bmy_rmg(:,:,1)));

end
