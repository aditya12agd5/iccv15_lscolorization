function fmat = pyr_to_feature_vec(pyr, mask)

count = 1;
if(nargin < 2)
  for i = 1:numel(pyr)
    image = pyr{i};
    for j = 1:size(image,3)
      fmat(:, count) = reshape(image(:,:,j), [], 1);
      count = count + 1;
    end
  end
else
  for i = 1:numel(pyr)
    image = pyr{i};
    for j = 1:size(image,3)
      wl = image(:,:,j);
      fmat(:, count) = reshape(wl(mask), [], 1);
      count = count + 1;
    end
  end
end

end
