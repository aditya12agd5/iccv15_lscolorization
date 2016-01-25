function [xOpt, muOpt, Sigma, wOpt, bDistInit] = optimizeBhattacharyyaDist(x, muS, SigmaS, wS, muT, SigmaT, wT, Amat, bVec, wmbv, pcell, filters, nums, numt, deps, nchannels)

if SigmaS ~= SigmaT
	disp(['Variances of source and target should be same!']);
end

startDist = evalBhattacharyyaDist(muS, muT, wS, wT, SigmaS);
disp(['Initial Bhattacharyya Distance: ' num2str(startDist)]);

bDistInit = single(startDist);

startLearchDist = evalLearchObjective(x, Amat, wmbv, pcell, filters, nums, numt, deps);
disp(['Initial Learch Objective: ' num2str(startLearchDist)]);

N = size(x, 2);
M = size(muS, 2);
Sigma = SigmaS;

%lambda = 0.5; %weight between learch and bhattacharyya objective
lambda = 5; %weight between learch and bhattacharyya objective
%lambda = 0; %don't consider learch at all

dists = [startDist];
k = 1;

derivative_norm = 10^5;

while k < 10 & derivative_norm > 10^-5
	
	%compute the gradient direction
	MtM_x = applymat_filter(Amat, pcell, filters, nums, numt, deps, reshape(x, 1, [])');		
	
	%Double derivative gradient update
	%MtM2_x = applymat_filter(Amat, pcell, filters, nums, numt, deps, MtM_x);
	%MtM_Mtb = applymat_filter(Amat, pcell, filters, nums, numt, deps, bVec);
	%deltaXLearch = reshape((MtM2_x-(2.*MtM_Mtb)), nchannels, []);

	deltaXLearch = reshape(((2.*MtM_x)-(2.*bVec)), nchannels, []);

	for j = [1:N]
		xCurr = x(:, j);
		n_ij = (1/(2*pi*Sigma)).*exp((-1/(2*Sigma^2)).*(sum((repmat(xCurr, [1 M])-muS).^2)));
		deltaX(:, j) = (0.25/(Sigma^2)).*sum(((muS-muT).*repmat(n_ij, [2 1])), 2);
	end

	deltaX = deltaX + lambda.*(deltaXLearch);

	%perform line search for step size
	alpha = lineSearch(x, deltaX, muS, muT, wS, wT, Sigma, Amat, wmbv, pcell, filters, nums, numt, deps, lambda);

	%update x
	x = x - alpha.*deltaX;

	derivative_norm = alpha.*sqrt(sum(sum((deltaX.*deltaX))))/N;

	%updated bhattacharyya distance, muS, wS
	[bDist, muS, wS] = evalBhattacharyyaDistX(x, muS, muT, wS, wT, Sigma);
	bLearch = evalLearchObjective(x, Amat, wmbv, pcell, filters, nums, numt, deps);
	
	disp(['Iteration# ' num2str(k) ' Derivative Norm: ' num2str(derivative_norm)]);
	disp(['Iteration# ' num2str(k) ' Bhattacharyya Distance: ' num2str(bDist) ' Learch Objective: ' num2str(bLearch)]);
	bDist = bDist + lambda.*bLearch;

	dists(end+1) = bDist;
	
	k = k + 1;
end

xOpt = x;
muOpt = muS;
wOpt = wS;

end

function alpha = lineSearch(x, deltaX, muS, muT, wS, wT, Sigma, Amat, wmbv, pcell, filters, nums, numt, deps, lambda)

a=0;                            % start of interval
b=0.1;                            % end of interval
epsilon=10^-5;               % accuracy value
iter= 12;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations


x1=tau*a+(1-tau)*b;             % computing x values
x2=(1-tau)*a+tau*b;

%call bhattacharyya dist
[f_x1, ~, ~] = evalBhattacharyyaDistX(x - x1.*deltaX, muS, muT, wS, wT, Sigma);                     
f_x1 = f_x1 + lambda*evalLearchObjective(x - x1.*deltaX, Amat, wmbv, pcell, filters, nums, numt, deps);

[f_x2, ~, ~] = evalBhattacharyyaDistX(x - x2.*deltaX, muS, muT, wS, wT, Sigma);
f_x2 = f_x2 + lambda*evalLearchObjective(x - x2.*deltaX, Amat, wmbv, pcell, filters, nums, numt, deps);


%disp('----------------- Start Line Search ---------------------');
while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;

    if(f_x1<f_x2)
        b=x2;
        x2=x1;
        x1=tau*a+(1-tau)*b;
        
	%call bhattcharyya dist
	f_x2 = f_x1;
	[f_x1, ~, ~] = evalBhattacharyyaDistX(x - x1.*deltaX, muS, muT, wS, wT, Sigma);                     
	f_x1 = f_x1 + lambda*evalLearchObjective(x - x1.*deltaX, Amat, wmbv, pcell, filters, nums, numt, deps);
        
    else
        a=x1;
        x1=x2;
        x2=(1-tau)*a+tau*b;
        
	%call bhattcharyya dist
	f_x1 = f_x2;
	[f_x2, ~, ~] = evalBhattacharyyaDistX(x - x2.*deltaX, muS, muT, wS, wT, Sigma);
	f_x2 = f_x2 + lambda*evalLearchObjective(x - x2.*deltaX, Amat, wmbv, pcell, filters, nums, numt, deps);
    end
    
end

% chooses minimum point
if(f_x1<f_x2)
%    sprintf('alpha_min=%f', x1)
%    sprintf('Bhattacharyya Dist=%f ', f_x1)
    alpha = x1;
else
%    sprintf('alpha_min=%f', x2)
%    sprintf('Bhattacharyya Dist=%f ', f_x2)
    alpha = x2;
end

end

function dist = evalLearchObjective(x, Amat, wmbv, pcell, filters, nums, numt, deps)
	%Derivative objective
	%MtM_x = applymat_filter(Amat, pcell, filters, nums, numt, deps, reshape(x, 1, [])');
	%dist = (1/(nums*numt))*sum((MtM_x - bVec).^2);

	vec = reshape(x, 1, [])';	
	im = reshape(reshape(vec, [], nums*numt)', nums, numt, []);
	filtered_im = applyfilters(im, pcell, filters);
	im_vec = reshape(reshape(filtered_im,[],size(filtered_im,3))',[],1);
	%size_amat = size(Amat)
	%size_imvec = size(im_vec)
	%size_wmbv = size(wmbv)
	dist = (1/size(wmbv,1))*sum((Amat*im_vec - wmbv).^2);
end


function [dist, muS, wS] = evalBhattacharyyaDistX(x, muS, muT, wS, wT, Sigma)

	N = size(x, 2);
	M = size(muS, 2);
	d = size(x, 1);

	wSUpdate = double(zeros(size(wS)));
	eps = 10^-9;

	% start opt
	n_ij = double(zeros(N, M));
	for i = [1:M]
		n_ij(:, i) =(wS(i)/(2*pi*Sigma)).*exp((-1/(2*Sigma^2)).*(sum((x - repmat(muS(:, i), [1 N])).^2, 1)')) ;
	end

	sum_n_ij = sum(n_ij, 2);
	%sum_n_ij(sum_n_ij < eps) = eps;
	%n_ij(n_ij < eps) = eps;

	n_ij = n_ij./repmat(sum_n_ij, [1 M]);
	n_ij(find(isnan(n_ij) == 1)) = 0;
	
	wSUpdate = sum(n_ij, 1);
	muS = (x*n_ij)./(repmat(wSUpdate, [d 1]));
	wS = wSUpdate./N;
	% end opt

	dist = evalBhattacharyyaDist(muS, muT, wS, wT, Sigma);
end

function dist = evalBhattacharyyaDist(muS, muT, wS, wT, Sigma)
	dist = ((0.125/(Sigma^2))*sum(sum((muS-muT).^2))) - (0.5*sum(log(wS.*wT)));
end


