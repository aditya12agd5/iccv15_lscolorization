function [x, h, xdata] = obtainbmyrmghist(I)

eps = 0.0001;

if(size(I,3) == 3) %RGB Image
	I_gray = (I(:,:,1) + I(:,:,2) + I(:,:,3))./3;
	I_bmy = (I(:,:,3)-(0.5.*(I(:,:,1)+I(:,:,2))))./(I_gray+eps);
	I_rmg = (I(:,:,1)-I(:,:,2))./(I_gray+eps);
else %BMY RMG Image
	I_bmy = I(:, :, 1);
	I_rmg = I(:, :, 2);
end

edges{1} = [-3:0.05:3];
edges{2} = [-3:0.05:3];

xdata = [I_bmy(:) I_rmg(:)]';
h = hist3([I_bmy(:) I_rmg(:)], 'Edges', edges);

x = [];
for i = [1:size(h,1)]
	for j = [1:size(h,2)]
		xCord = -3+((i-1)*0.05);
		yCord = -3+((j-1)*0.05);
		if(h(i,j) > 0)
			x(end+1:end+h(i,j), :) = repmat([xCord yCord], h(i,j), 1);
		end
	end
end

end
