function [ inDiv, outDiv ] = diversity( wtMat, PrefFRs, show )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [M, N] = size(wtMat);
    wtMat = abs(wtMat); %No negative elements allowed
    [kIn, kOut, ~] = nodeDegrees(wtMat);
    
    inNorm = wtMat ./ repmat(sum(wtMat), M, 1);
    inNorm = inNorm .* full(spfun(@(x) log2(x), inNorm));
    inDiv = -sum(inNorm) ./ log2(kIn);
    
    outNorm = wtMat ./ repmat(sum(wtMat, 2), 1, N);
    outNorm = outNorm .* full(spfun(@(x) log2(x), outNorm));
    outDiv = -sum(outNorm, 2)' ./ log2(kOut);

    if show == 1
       figure; hist(inDiv, N/10 + 1); title('In Diversity');
       figure; hist(outDiv, M/10 + 1); title('Out Diversity');
       figure; scatter(PrefFRs, inDiv); title('PFRs v. In Div');
       figure; scatter(PrefFRs, outDiv); title('PFRs v. Out Div');  
    end
    
end

