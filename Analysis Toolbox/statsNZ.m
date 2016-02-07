function [ sumNZ, meanNZ, stdNZU, stdNZL ] = statsNZ( wtMat, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [~, el] = size(wtMat);
    sumNZ = zeros(1, el);
    meanNZ = zeros(1, el);
    stdNZU = zeros(1, el);
    stdNZL = zeros(1, el);
    
    for i = 1:el
        if isempty(varargin)
            col = nonzeros(wtMat(:, i));
        else
            col = wtMat(varargin{1}(:, i), i);
        end
       if numel(col) ~= 0
           sumNZ(1, i) = sum(col);
           meanNZ(1, i) = mean(col);
           [stdNZL(1, i), stdNZU(1, i)]= semistd(col);
       end
        
    end

end

