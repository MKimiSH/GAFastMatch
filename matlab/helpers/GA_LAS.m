function [oinds, dist] = GA_LAS(inds, dist, sigma)

mindist = min(dist);
maxdist = max(dist);

% note that _inds_ may not be successive
% and _indpool_ corresponds to indices of _dist_ 
% so samples(inds(indpool==i)) corresponds to samples in this bin
% rather than samples(indpool==i)
thres = linspace(mindist, maxdist*(1.01), sigma+1);
indpool = zeros(length(inds),1);
ninds = zeros(sigma,1);
for i = 1:sigma
    curinds = find(dist<thres(i+1) & dist>=thres(i));
    ninds(i) = length(curinds);
    indpool(curinds) = i;
end

% nlvl = min(ninds);
nlvl = round(length(inds)/sigma);
oinds = [];

for i = 1:sigma
    curinds = find(indpool==i);
    %cursamples = curinds(randperm(ninds(i), nlvl));
    if(~isempty(curinds))
        cursamples = curinds(randi(ninds(i), nlvl, 1));
        oinds = [oinds, inds(cursamples)];
    end
    
end

oinds = oinds(randperm(length(oinds)));

end