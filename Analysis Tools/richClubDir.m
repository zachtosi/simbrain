function [ rccs, richness, significance ] = richClubDir( adjMat, riParam )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [n, ~] = size(adjMat);
    maxVal = max(riParam);
    [rccs, richness] = rc_dir(adjMat, riParam);    
    rccShuff = zeros(20, numel(1:maxVal/200:maxVal));
    richnessShuff = zeros(20, n);
    for i = 1:20
        [rccShuff(i, :), richnessShuff(i, :)] = rc_dir(dir_generate_srand(adjMat), riParam);
    end
    plot(1:maxVal/200:maxVal, rccs, 'r');   
    meanRccSh = mean(rccShuff);
    meanRichnessSh = mean(richnessShuff);
    stdRccSh = std(rccShuff);
    hold; errorbar(1:maxVal/20:maxVal, ...
        meanRccSh(1, 1:10:200), ...
         stdRccSh(1, 1:10:200));
     
    plot(1:maxVal/200:maxVal, rccs ./ meanRccSh, 'k');
    mx = 1:maxVal/20:maxVal;
    plot(mx, ones(1, numel(mx)), '--k');
    hold;
    significance = abs(richness - meanRichnessSh) ./ std(richnessShuff);
    richness = richness ./ meanRichnessSh;
    rccs = rccs ./ meanRccSh;
end

