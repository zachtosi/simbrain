function [ out, outRw ] = neighborHist ( wts )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    binMat = wts ~= 0;
    rw = dir_generate_srand(binMat);
    binMat = (binMat + binMat') ~= 0;
    rw = (rw + rw') ~= 0;
    [n, m] = size(wts);
    out = zeros(n, m);
    outRw = zeros(n, m);
    for i = 1:n
        neigh = find(binMat(i, :) ~= 0);
        neighRw = find(rw(i, :) ~= 0); 
        nn = numel(neigh);
        nnRw = numel(neighRw);
        for j = 1:nn
            k = neigh(j);
            if (k > i)
                out(i, k) = sum(binMat(i, :) .* binMat(k, :));    
            end
        end
        for j = 1:nnRw
            k = neighRw(j);
            if (k > i)
                outRw(i, k) = sum(rw(i, :) .* rw(k, :));    
            end
        end
    end
    
    ncs = nonzeros(out)';
    ncsRw = nonzeros(outRw)';
    
    m1 = max(ncs);
    m2 = max(ncsRw);
    if m1 > m2
        MAX = m1;
        figure; hist(ncs, MAX/2); title('Model');
        figure; hist([ ncsRw, MAX ], MAX/2); title('Null');
    else
        MAX = m2;
        figure; hist([ ncs, MAX ], MAX/2); title('Model');
        figure; hist(ncsRw, MAX/2); title('Null'); 
    end

    

end

