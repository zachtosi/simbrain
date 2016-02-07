function [ rccs, richness, significance ] = richClubDir( adjMat, riParam )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [n, ~] = size(adjMat);
    maxVal = max(riParam);
    [rccs, richness] = rc_dir(adjMat, riParam);    
    rccShuff = zeros(100,n);
    richnessShuff = zeros(100, n);
    parfor i = 1:100
        [rccShuff(i, :), richnessShuff(i, :)] = rc_dir(dir_generate_srand(adjMat), riParam);
    end
    sortedRIP = sort(riParam, 'ascend');
    plot(sortedRIP, rccs, 'r');   
    meanRccSh = mean(rccShuff);
    meanRichnessSh = mean(richnessShuff);
   [stdRccShL, stdRccShU] = semistd(rccShuff);
     hold; shadedErrorBar(sortedRIP, ...
        meanRccSh, ...
         [stdRccShU; stdRccShL], 'b', 1);
     
    plot(sortedRIP, rccs ./ meanRccSh, 'k');
    mx = 1:maxVal/20:maxVal;
    plot(mx, ones(1, numel(mx)), '--k');
    hold;
    significance = abs(richness - meanRichnessSh) ./ std(richnessShuff);
    richness = riParam .* ((richness ./ meanRichnessSh));
    rccs = rccs ./ meanRccSh;
end

