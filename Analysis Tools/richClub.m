function [rccs] = richClub( adjMat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    adjMat = adjMat ~= 0;
    degs = sum(adjMat);
    maxDeg = max(degs);
    
    rccs = richclubcoeff(adjMat, 1:maxDeg);
    rccShuff = zeros(20, maxDeg);
    for i = 1:20
        rccShuff(i, :) = richclubcoeff(dir_generate_srand(adjMat),...
            1:maxDeg);
    end
    plot(rccs, 'r');   
    meanRccSh = mean(rccShuff);
    stdRccSh = std(rccShuff);
    hold; errorbar([1 int32(ceil(maxDeg/20)):int32(ceil(maxDeg/20)):maxDeg maxDeg], ...
        meanRccSh(1, [1 int32(ceil(maxDeg/20)):int32(ceil(maxDeg/20)):maxDeg maxDeg]), ...
         stdRccSh(1, [1 int32(ceil(maxDeg/20)):int32(ceil(maxDeg/20)):maxDeg maxDeg]));
     
    plot(rccs ./ meanRccSh, 'k');
    mx = [1:int32(ceil(maxDeg/20)):maxDeg maxDeg];
    plot(mx, ones(1, numel(mx)), '--k');hold;
end

