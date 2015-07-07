function [ vers, degs ]...
    = nodeVersatility( wtMat )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    ei = sum(wtMat, 2)' > 0;
    len = length(wtMat);
    vers = zeros(1, len);
    binaryMat = wtMat ~= 0;
    degs = sum(binaryMat, 1) + sum(binaryMat, 2)';
    
    for i = 1:len
       [~, out, ~] = find(binaryMat(i, :) ~= 0);
       [in, ~, ~] = find(binaryMat(:, i) ~= 0);
       neighborDegs = zeros(1, length(in) + length(out));
       for p = 1:length(in)
          neighborDegs(1, p) = degs(1, in(p)); 
       end
       for q = 1:length(out)
          neighborDegs(1, q+length(in)) = degs(1, out(q)); 
       end
       neighborDegs = abs(neighborDegs - degs(i));
       neighborDegs = neighborDegs - mean(neighborDegs);
       neighborDegs = neighborDegs .* neighborDegs;
       vers(i) = sqrt(sum(neighborDegs)/degs(i));
    end
    figure;
    scatter(degs, vers);
    %exDegs = nonzeros(degs .* ei);
    %exVers = nonzeros(vers .* ei);
    %scatter(exDegs, exVers, 'r');
    %hold;
    %inDegs = nonzeros(degs .* ~ei);
    %inVers = nonzeros(vers .* ~ei);
    %scatter(inDegs, inVers, 'b');
end

