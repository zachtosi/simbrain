function [ sumNZ, meanNZ, stdNZ ] = statsNZ( wtMat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [~, el] = size(wtMat);
    sumNZ = zeros(1, el);
    meanNZ = zeros(1, el);
    stdNZ = zeros(1, el);
    
    for i = 1:el
       col = nonzeros(wtMat(:, i));
       if numel(col) ~= 0
           sumNZ(1, i) = sum(col);
           meanNZ(1, i) = mean(col);
           stdNZ(1, i) = std(col);
       end
    end

end

