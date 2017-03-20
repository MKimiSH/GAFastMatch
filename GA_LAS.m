function [oinds, dist] = GA_LAS(inds, dist, sigma)

mindist = min(dist);
maxdist = max(dist);

% note that _inds_ may not be successive
% and _indpool_ corresponds to indices of _dist_ 
% so samples(inds(indpool==i)) corresponds to samples in this bin
% rather than samples(indpool==i)
thres = linspace(mindist, maxdist+0.1, sigma+1);
indpool = zeros(size(inds),1);
ninds = zeros(sigma,1);
for i = 1:sigma
    curinds = find(dist<thres(i+1) & dist>=thres(i));
    ninds(i) = length(curinds);
    indpool(curinds) = i;
end

nlvl = min(ninds);
oinds = [];

for i = 1:sigma
    curinds = find(indpool==i);
    cursamples = curinds(randperm(ninds, nlvl));
    oinds = [oinds; inds(cursamples)];
end

end