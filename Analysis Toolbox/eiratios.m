function [ p, ee, ei, ie, ii ] = eiratios( wtMat )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    [m, n] = size(wtMat);
    bin = wtMat ~= 0;
    p = sum(sum(bin)) / (m * (n-1));
    exin = sum(wtMat, 2) > 0;

    ee = sum(sum(bin(exin, exin))) / (sum(exin)*(sum(exin)-1));
    ei = sum(sum(bin(exin, ~exin))) / (sum(exin) * sum(~exin));
    ie = sum(sum(bin(~exin, exin))) / (sum(exin) * sum(~exin));
    ii = sum(sum(bin(~exin, ~exin))) / (sum(~exin) * (sum(~exin)-1));

end

