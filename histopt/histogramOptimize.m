function [res, corresp, bDistInit] = histogramOptimize(estImgVec, Amat, bvec, wmbv, pcell, filters, width, height, deps, exemplarHist, nchannels, bandwidth, image_prefix)

addpath(genpath('./gmm_code/'));

x = histToData(exemplarHist); 

[clustCent, point2cluster, clustMembsCell] = MeanShiftClusterGaussian(x', bandwidth); 
numClust = length(clustMembsCell);
[gtClustCent, selectMembers] = thresholdClusters(x', clustCent, clustMembsCell, numClust, 1, image_prefix);

em_gmm_params = initEmGmmParams(gtClustCent, bandwidth/2);
gmmGt = EM_GM_fast(x, size(gtClustCent, 2), [], [], 0, em_gmm_params);
 	
muTarget = gmmGt.mu';
SigmaTarget = bandwidth/2;
wTarget = gmmGt.PComponents;

estImg = reshape(reshape(estImgVec, nchannels, [])', width, height, []);
[x, h, xdata] = obtainbmyrmghist(estImg);


[clustCent, point2cluster, clustMembsCell] = MeanShiftClusterGaussian(xdata, bandwidth, gtClustCent); 
numClust = length(clustMembsCell);
disp(['#clusters after meanshift on est: ' num2str(numClust)]);
[~, selectMembers, validCorrespClusters] = thresholdClusters(xdata, clustCent, clustMembsCell, numClust, 2, image_prefix);

disp(['#clusters that found correspondence: ' num2str(size(validCorrespClusters, 2))]);

muTarget = muTarget(:, validCorrespClusters); 
wTarget = wTarget(validCorrespClusters); 

wTarget = wTarget./sum(wTarget);

clustCent = clustCent(:, validCorrespClusters); 

em_gmm_params = initEmGmmParams(clustCent, bandwidth/2);
gmmEst = EM_GM_fast(xdata(:, selectMembers)', size(clustCent, 2), [], [], 0, em_gmm_params);
	

muSource = gmmEst.mu';
SigmaSource = bandwidth/2;
wSource = gmmEst.PComponents; 

[xOptimizedData, muOpt, SigmaOpt, wOpt, bDistInit] = optimizeBhattacharyyaDist(xdata, muSource, SigmaSource, wSource, ...
	muTarget, SigmaTarget, wTarget, Amat, bvec, wmbv, pcell, filters, width, height, deps, nchannels);

res = reshape(xOptimizedData, 1, [])';
corresp = numel(validCorrespClusters);

if(size(res) ~= size(estImgVec))
	disp(['ERROR in return vec size']);
end

end

function [newClustCent, selectMembers, validIds] = thresholdClusters(x, clustCent, clustMembsCell, numClust, id, image_prefix)

cVec = 'bgrcmybgrcmybgrcmybgrcmy'; cVec = [cVec cVec];
countcVec = 1;

newClustCent = [];
selectMembers = [];

validIds = [];
for k = 1:min(numClust,length(cVec))
    myMembers = clustMembsCell{k};
    if numel(myMembers) > 1000 | id ~= 1
	    myClustCen = clustCent(:,k);
	    newClustCent = horzcat(newClustCent, myClustCen);  
	    if(numel(myMembers) > 25)
		  validIds(end+1) = k;
	    end
	    countcVec = countcVec + 1;
	    selectMembers = [selectMembers myMembers];
   end
end

end

function em_gmm_params = initEmGmmParams(clustCent, bandwidth)

em_gmm_params = {};
d = size(clustCent, 1);
k = size(clustCent, 2); 
em_gmm_params.M = clustCent;
em_gmm_params.W = ones(1, k)./k;
em_gmm_params.V = double(zeros(d, d, k));
for nGmm=[1:k]
	em_gmm_params.V(:, :, nGmm) = (bandwidth^2).*eye(d, d);
end

end

function [] = plotGmm(gmmObj, id, image_prefix)

f2 = figure(3);

clf, hold on;
for itrGMMs = [1:size(gmmObj.mu, 1)]
	currObj = gmdistribution(gmmObj.mu(itrGMMs, :), gmmObj.Sigma(:, :, itrGMMs), 1);
	ezcontour(@(p,q)pdf(currObj,[p q]),[-3 3],[-3 3], 200);
end
hold off;
title(['Weights ' num2str(gmmObj.PComponents(:)')]);
saveas(f2, [image_prefix '_gmm-' num2str(id) '.png']);

end
