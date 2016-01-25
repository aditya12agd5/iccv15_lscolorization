function flip_kernel = doubleflip(kernel)
assert(size(kernel,1) == size(kernel,2));

flip_kernel = zeros(size(kernel));
for i = 1:size(kernel,1)
  for j = 1:size(kernel,2)
    flip_kernel(i,j) = kernel(end-i+1, end-j+1);
  end
end
end
