function [clustCent,data2cluster,cluster2dataCell] = MeanShiftCluster(dataPts,bandWidth, initMeans, plotFlag);
%(@aditya) Modified to use a Gaussian Kernel
%perform MeanShift Clustering of data using a flat kernel
%
% ---INPUT---
% dataPts           - input data, (numDim x numPts)
% bandWidth         - is bandwidth parameter (scalar)
% plotFlag          - display output if 2 or 3 D    (logical)
% ---OUTPUT---
% clustCent         - is locations of cluster centers (numDim x numClust)
% data2cluster      - for every data point which cluster it belongs to (numPts)
% cluster2dataCell  - for every cluster which points are in it (numClust)
% 
% Bryan Feldman 02/24/06
% MeanShift first appears in
% K. Funkunaga and L.D. Hosteler, "The Estimation of the Gradient of a
% Density Function, with Applications in Pattern Recognition"


%*** Check input ****
if nargin < 2
    error('no bandwidth specified')
end

if nargin < 4
    plotFlag = false;
else
    plotFlag = true;
end

%**** Initialize stuff ***
[numDim,numPts] = size(dataPts);
numClust        = 0;
bandSq          = bandWidth^2;
initPtInds      = 1:numPts;
maxPos          = max(dataPts,[],2);                          %biggest size in each dimension
minPos          = min(dataPts,[],2);                          %smallest size in each dimension
boundBox        = maxPos-minPos;                        %bounding box size
sizeSpace       = norm(boundBox);                       %indicator of size of data space
stopThresh      = 1e-3*bandWidth;                       %when mean has converged
clustCent       = [];                                   %center of clust
beenVisitedFlag = zeros(1,numPts,'uint8');              %track if a points been seen already
numInitPts      = numPts;                               %number of points to posibaly use as initilization points
clusterVotes    = zeros(1,numPts,'double');             %used to resolve conflicts on cluster membership

countInitMeans = 1;


countInitMeans = 1;

if nargin < 3
numInitMeans = numPts;
else
disp(['Init Means# ' num2str(size(initMeans, 2))]);
numInitMeans = size(initMeans, 2);
end

while numInitPts & countInitMeans <= numInitMeans

    tempInd         = ceil( (numInitPts-1e-6)*rand);        %pick a random seed point

    if nargin < 3
	    stInd           = initPtInds(tempInd);                  %use this point as start of mean
	    myMean          = dataPts(:,stInd);                           % intilize mean to this points location
	    candMean = 0;
    else
	    %myMean = initMeans(:, countInitMeans);
	    candMean = initMeans(:, countInitMeans);
            [nearDist, nearInd] = min(sum((repmat(candMean,1,numPts) - dataPts).^2));    %dist squared from mean to all points still active
	    myMean = dataPts(:, nearInd);
	    countInitMeans  = countInitMeans + 1;
    end

    myMembers       = [];                                   % points that will get added to this cluster                          
    thisClusterVotes    = zeros(1,numPts,'double');         %used to resolve conflicts on cluster membership

    while 1     %loop untill convergence
        
        sqDistToAll = sum((repmat(myMean,1,numPts) - dataPts).^2);    %dist squared from mean to all points still active
        inInds      = find(sqDistToAll < bandSq);               %points within bandWidth
	
	if isnan(myMean)
		disp('!!!!!!!!!!!! Error !!!!!!!!!!!!!!!!!');
	end

	
	%@aditya: changed to smooth voting
	%thisClusterVotes(inInds) = thisClusterVotes(inInds)+1;  %add a vote for all the in points belonging to this cluster
       	weights = gaussianWeight(dataPts(:, inInds), myMean, bandWidth);	
	thisClusterVotes(inInds) = thisClusterVotes(inInds) + weights;
 
        
        myOldMean   = myMean;                                   %save the old mean
        %@aditya: changed to smooth voting
	%myMean      = mean(dataPts(:,inInds),2);                %compute the new mean
	myMean = sum(dataPts(:, inInds).*repmat(weights, [size(dataPts, 1) 1]), 2)/sum(weights);		

        myMembers   = [myMembers inInds];                       %add any point within bandWidth to the cluster
        beenVisitedFlag(myMembers) = 1;                         %mark that these points have been visited
        
        %*** plot stuff ****
        if plotFlag
            figure(12345),clf,hold on
            if numDim == 2
                plot(dataPts(1,:),dataPts(2,:),'.')
                plot(dataPts(1,myMembers),dataPts(2,myMembers),'ys')
                plot(myMean(1),myMean(2),'go')
                plot(myOldMean(1),myOldMean(2),'rd')
                pause
            end
        end
        
	%**** if mean doesn't move much stop this cluster ***
        if norm(myMean-myOldMean) < stopThresh  | (candMean ~= 0 & norm(candMean-myMean) > bandWidth/2)

	    if (candMean ~= 0 & norm(candMean-myMean) > bandWidth/2)
            	%myMean = (candMean+myMean)/2;
		disp(['Stopped mode from moving too far from init']);
	    end

	    %check for merge posibilities
            mergeWith = 0;

	    if nargin < 3
  	          for cN = 1:numClust
  	              distToOther = norm(myMean-clustCent(:,cN));     %distance from posible new clust max to old clust max
  	              if distToOther < bandWidth/2                    %if its within bandwidth/2 merge new and old
  	                  mergeWith = cN;
 	                   break;
 	               end
 	           end
	    end
            
            
            if mergeWith > 0    % something to merge
                clustCent(:,mergeWith)       = 0.5*(myMean+clustCent(:,mergeWith));             %record the max as the mean of the two merged (I know biased twoards new ones)
                clusterVotes(mergeWith,:)    = clusterVotes(mergeWith,:) + thisClusterVotes;    %add these votes to the merged cluster
            else    %its a new cluster
                numClust                    = numClust+1;                   %increment clusters
                clustCent(:,numClust)       = myMean;                       %record the mean  
                clusterVotes(numClust,:)    = thisClusterVotes;
            end

            break;
        end

    end
    
    
    initPtInds      = find(beenVisitedFlag == 0);           %we can initialize with any of the points not yet visited
    numInitPts      = length(initPtInds);                   %number of active points in set
    if nargin >= 3
	numInitPts = numPts; 
    end
end

falseRow = zeros(1, size(clusterVotes, 2));
clusterVotes = vertcat(falseRow, clusterVotes);
[val,data2cluster] = max(clusterVotes,[],1);                %a point belongs to the cluster with the most votes
data2cluster = data2cluster - 1;

%*** If they want the cluster2data cell find it for them
if nargout > 2
    cluster2dataCell = cell(numClust,1);
    for cN = 1:numClust
        myMembers = find(data2cluster == cN);
        cluster2dataCell{cN} = myMembers;
    end
end

end

function weights = gaussianWeight(X, mu, sigma)
	sqDist = sum((X - repmat(mu, [1 size(X, 2)])).^2);
	weights = (1/(2*pi*sigma)).*exp((-0.5/(sigma^2)).*sqDist);
end
