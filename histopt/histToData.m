function [x] = histToData(h)

%flatten histogram to data
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
