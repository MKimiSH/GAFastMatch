function [inds, dist] = GA_ShrinkGroupbyLambda(inds, dist, lambda)

n = length(dist);
%% sample a subset to calculate threshold when necessary
pct = 0.1;
if n>1000 
    thresh = sampleThresh(dist, pct, lambda);
else
    sorteddist = sort(dist, 'ascend');
    thresh = sorteddist(round(n*lambda));
end
%% use thresh to choose samples
inliers = find(dist<thresh);
inds = inds(inliers);
dist = dist(inliers);

%% use the threshold on all distances
end

function [thresh] = sampleThresh(dist, pct, lambda)
n = length(dist);
nsample = round(n*pct);
indsample = randperm(n, nsample);
sample = dist(indsample);
sample = sort(sample, 'ascend');
thresh = sample(round(nsample*lambda));
end