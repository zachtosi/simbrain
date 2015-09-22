function [ rcc, richness ] = rc_dir( mat, riParam )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [n, ~] = size(mat);
    binMat = mat ~= 0;
    ma = max(riParam);
    rcc = zeros(1, numel(1:ma/200:ma));
    j = 1;
    richness = zeros(1, n);
    for i = 1:ma/200:ma
        gr = riParam > i;
        num = sum(gr);
        if num < 2
            break;
        end
        rcc(j) = sum(sum(binMat(gr, gr))) / (num * (num - 1));
        richness(gr) = rcc(j);
        j = j + 1;
    end
    
end

