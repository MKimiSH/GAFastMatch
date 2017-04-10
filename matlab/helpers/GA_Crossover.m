function [newsamples] = GA_Crossover(samples, n)

% nsamples = size(samples, 1);
% nbits = 6 * 2^n;
pcross = 0.7;

% randomly permute the samples
nsamples = size(samples,1);
order = randperm(nsamples);
samples = samples(order,:);
samples = samples(order,:);
samples = samples(order,:);

logicalSamples = GA_LogicalSamplesFromInt(samples, n);%zeros(nsamples, nbits);

logicalSamples = GA_LogicalCrossover(logicalSamples, pcross);

newsamples = GA_IntSamplesFromLogical(logicalSamples, n);


end

function [lsam] = GA_LogicalSamplesFromInt(samples, n)

nsamples = size(samples, 1);
nbits = 6 * n;

lsam = zeros(nsamples, nbits, 'logical');
for i=n:-1:1
    % n=3: lsam(1,:) = tx2,tx1,tx0,ty2,ty1,ty0,...,\the22,\the21,\the20
    lsam(:, i:n:nbits) = mod(samples, 2);
    samples = samples / 2;
end
end

function [lsam] = GA_LogicalCrossover(lsam, pcross)

[nsamples,nbits] = size(lsam);

for i=1:2:nsamples-1
    if(rand(1)<pcross)
        p1 = lsam(i, :);
        p2 = lsam(i+1, :);
        ptcross = randi(nbits);
        t = p1(ptcross:nbits);
        p1(ptcross:nbits) = p2(ptcross:nbits);
        p2(ptcross:nbits) = t;
        lsam(i, :) = p1;
        lsam(i+1,:) = p2;
    end
    
end
end

function [sam] = GA_IntSamplesFromLogical(lsam, n)

[nsamples,nbits] = size(lsam);
assert(nbits == n*6);
sam = zeros(nsamples, 6, 'uint8');

for i=1:n
    sam = sam * 2;
    sam = sam + uint8(lsam(:, i:n:nbits));
end

end