function [bestConfig,bestTransMat,sampledError] = ...
        GA_FindBestTransformation(I1,I2,bounds,steps,epsilon,delta,photometricInvariance,templateMask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% verify input image types
if ( ~strcmp(class(I1),'double') || ~strcmp(class(I2),'double')) %#ok<STISA>
    error('GAFastMatch: I1 and I2 should both be of class ''double'' (in the range [0,1])');
end

if ((size(templateMask,1) ~= size(I1,1)) || (size(templateMask,2) ~= size(I1,2)) )
    error('GAFastMatch: Template mask not same size as template');
end

% if (size(I1,1)~=size(I1,2) || size(I2,1)~=size(I2,2))
%     error('GAFastMatch: Template and image should be square!');
% end

isGrayscale = (size(I1,3)==1);

if ~isGrayscale
%     error('FastMatch: Only grayscale images are currently supported');
end


%% blur in main loop - this reduces the total-variation and gives better results
if (isGrayscale)
    origI1 = I1;
    origI2 = I2;
end
[h1,w1,d] = size(I1);

%% generate Theta(1/eps^2) random points (and fresh ones each iteration later on)
numPoints = round(10/epsilon^2);
[xs, ys] = getPixelSample(templateMask, numPoints);

%% generate the Net
% [configs,gridSize] = CreateListOfConfigs(bounds,steps);
% the command above will fail because of too many configs!!!!
% only create list of configs of those in the group!

% if (size(configs,1) > 71000000)
%         error('more than 35 million configs!');
% end


%% main loop
% TO Change!

% parameters
% n1 = n1; % template size 
% n2 = n2; % target image size, if the images are not square, reshape them! 
[h2,w2,d] = size(I2);
n = 8; % length of the code --> number of steps = 2^n
sigma = 5; % LAS sampling parameter 
eps = 3; % step for SAD computation
delta = 11^6; % initial group size
lambda = 0.9; % reduction of group size per generation
alpha = 0; % two parameters for to bound the lambda 
beta = 0; % of which I question the use
c = 200; % #samples of last generation

% procedure
% 1. initialization
% 2. for every generation until generation size is smaller than c
% <1. sample learning set
% <2. tune alpha and beta with learning set
% <3. use LAS to sample the next generation and reduce the size with lambda
% <4. crossover the set in <3.
% 3. return the best A in this small final group

% 1-tx, 2-ty, 3-theta2, 4-sx, 5-sy, 6-theta1; -temporary
% Dr. Zhang: 1-tx, 2-ty, 3-theta2, 4-theta1, 5-sx, 6-sy;
initSamples = randi([0, 2^n - 1], delta, 6, 'uint8');
% initConfigs = GA_CreateListOfConfigs(bounds,steps,initSamples);
group = ones(delta,1);
groupSize = delta;
samples = initSamples;
% configs = initConfigs;
niter = 0;
winnerAffines = [];
championAffine = [];
totTime = 0;
while(groupSize > c)
    niter = niter + 1;
    % only evaluate SAD of configs ONCE
    fprintf('iteration %d\n', niter);
    % oldsamples -> configs
    configs = GA_CreateListOfConfigs(bounds,steps,samples);
    groupSize = size(samples, 1); % #samples this iteration
    groupidx = 1:groupSize;
    r1x = 0.5*(w1-1);
    r1y = 0.5*(h1-1);
    r2x = 0.5*(w2-1);
    r2y = 0.5*(h2-1);
%     groupConfigs = configs(groupidx);
    
    % 2] configs -> affine exclude outliers
    Configs2AffineMEX = tic;
    [matrixConfigs_mex, insiders] = ...
        Configs2Affine_mex(configs',int32(h1), int32(w1), int32(h2), int32(w2), int32(r1x), int32(r1y), int32(r2x), int32(r2y));
    inBoundaryInds = find(insiders);
    matrixConfigs_mex = matrixConfigs_mex(:,inBoundaryInds);
    origNumConfigs = size(configs,1);
    configs = configs(inBoundaryInds,:);
    groupidx = groupidx(inBoundaryInds);
    samples = samples(inBoundaryInds,:);
    groupSize = length(groupidx);
    groupidx = 1:groupSize;
    Configs2Affine_mex_time = toc(Configs2AffineMEX);
    
    % 3] affine -> distances
    % distances is a groupSize x 1 vector
    EvaluateConfigsMEX = tic;
    if (isGrayscale)
        distances = EvaluateConfigs_mex(I1',I2',matrixConfigs_mex,int32(xs),int32(ys),int32(photometricInvariance));
        fprintf('----- GA: Evaluate Configs grayscale, with %d configs -----\n',size(configs,1));
    else
        distances = EvaluateConfigsVectorized_mex(permute(I1,[3,2,1]),permute(I2,[3,2,1]),matrixConfigs_mex,int32(xs),int32(ys),int32(photometricInvariance));
        fprintf('----- GA: Evaluate Configs vectorized, with %d configs -----\n',size(configs,1));
    end
    EvaluateConfigs_mex_time = toc(EvaluateConfigsMEX);
    totTime = totTime + Configs2Affine_mex_time + EvaluateConfigs_mex_time;
    
    % 4] use distances to shrink group with lambda
    % groupidx and distances are shrinked
    [groupidx, distances] = GA_ShrinkGroupbyLambda(groupidx, distances, lambda);

    % 5] use remaining distances for LAS
    [groupidx, distances] = GA_LAS(groupidx, distances, sigma);
    groupSize = length(groupidx);
    % groupidx -> newsamples
    samples = samples(groupidx, :);
    fprintf('mindist = %.3f\n', min(distances));
    
    % if crossover is allowed then optimality will be severely damaged.
    % 6] if small group then choose best one and stop.
    if(groupSize <= c)
        winnerAffines = matrixConfigs_mex(:, groupidx);
        [bestDist, bdidx] = min(distances);
        bestConfig = configs(bdidx, :);
        bestTransMat = CreateAffineTransformation(bestConfig); 
%         [winners, champion] = ...
%             GA_FindChampion(samples, distances, configs, matrixConfigs_mex);
        fprintf('stopping!!\n');
        fprintf('$$$ bestDist = %.3f\n', bestDist);
        break;
    end
    
    % 7] newsamples --crossover--> next iteration.
    [samples] = GA_Crossover(samples, n);
    [xs, ys] = getPixelSample(templateMask, numPoints);
end
% bestConfig = GA_FindBestConfig(I1, I2, configs);


%% debug error
if isGrayscale
    [xs,ys] = meshgrid(1:w1,1:h1);
    xs = xs(:)';
    ys = ys(:)';
    BD_full = EvaluateConfigs_mex(I1',I2',bestTransMat([1 4 7 2 5 8])',int32(xs),int32(ys),int32(photometricInvariance));
    BD_full_orig = EvaluateConfigs_mex(origI1',origI2',bestTransMat([1 4 7 2 5 8])',int32(xs),int32(ys),int32(photometricInvariance));
    fprintf('bestDist: %.4f, BD_full: %.4f, BD_full_orig: %.4f\n',bestDist,BD_full,BD_full_orig);
else
    [xs,ys] = meshgrid(1:w1,1:h1);
    xs = xs(:)';
    ys = ys(:)';
    BD_full = EvaluateConfigsVectorized_mex(permute(I1,[3,2,1]), permute(I2,[3,2,1]), bestTransMat([1 4 7 2 5 8])',int32(xs),int32(ys),int32(photometricInvariance));
    % BD_full_orig = EvaluateConfigsVectorized_mex(origI1',origI2',bestTransMat([1 4 7 2 5 8])',int32(xs),int32(ys),int32(photometricInvariance));
    fprintf('Color image: bestDist: %.4f, BD_full: %.4f\n',bestDist,BD_full);
end

%% for output
sampledError = bestDist;

return



function [res,i] = IsMemberApprox(A,row,err)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = 0;
for i = 1 : size(A,1)
        if (norm(A(i,:)-row) < err)
                res = 1;
                return
        end
end


function [goodConfigs,tooHighPercentage,extremelyHighPercentage,veryLowPercentage,orig_percentage,thresh] = ...
    GetGoodConfigsByDistance(configs,bestDist,newDelta,distances,bestGridVec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% targetNum = 20000;
% thresh = bestDist + newDelta/3;
thresh = bestDist + GetThreshPerDelta(newDelta);
goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
numGoodConfigs = size(goodConfigs,1);
orig_percentage = numGoodConfigs/size(configs,1);

% too many good configs - reducing threshold
while (numGoodConfigs > 27000)
        thresh = thresh * 0.99;
        goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
        numGoodConfigs = size(goodConfigs,1);
end

if (isempty(goodConfigs))
         thresh = min(distances);
        goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
        if (size(goodConfigs,1)>10000)
                inds = find(distances <= thresh);
                goodConfigs = configs(inds(1:100), :); % all with the same error exactly - probably equivalent
        end
end  

tooHighPercentage = (orig_percentage > 0.05);
veryLowPercentage = (orig_percentage < 0.01);
extremelyHighPercentage = (orig_percentage > 0.2);

if (~isempty(bestGridVec))
    [exists,bestGridInd] = IsMemberApprox(goodConfigs,bestGridVec,1000*eps);
    if (~exists)
        disp('problem with configs');
    end
end


function [xs, ys] = getPixelSample(mask, numPoints)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
locs = find(mask);
ind = randi(length(locs), [1,numPoints]);
[ys,xs] = ind2sub(size(mask),locs(ind));
ys = ys';
xs = xs';


