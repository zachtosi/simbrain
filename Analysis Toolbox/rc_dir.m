function [ rcc, richness ] = rc_dir( mat, riParam )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [n, ~] = size(mat);
    binMat = mat ~= 0;
    rcc = zeros(1, n);
    j = 1;
    richness = zeros(1, n);
    sortedRIP = sort(riParam, 'ascend');
    for i = 1:numel(riParam)
        gr = riParam >= sortedRIP(i);
        num = sum(gr);
        if num < 2
            rcc(j) = 1;
            continue;
        end
        rcc(j) = sum(sum(binMat(gr, gr))) / (num * (num - 1));
        richness(gr) = rcc(j);
        j = j + 1;
    end
    
end

