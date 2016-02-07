function [avgOIM, stdOIM] = sigScatterInVOut( wts, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[~, n] = size(wts);
ei = (sum(wts, 2) > 0)';
error(nargchk(1, 4, nargin));
sigP = 0.1;
if length(varargin) == 0
   NUM_NULL = 100;
elseif length(varargin) == 1
    NUM_NULL = varargin{1};
elseif length(varargin) == 2
    NUM_NULL = varargin{1};
    sigP = varargin{2};
elseif length(varargin) == 3
    NUM_NULL = varargin{1};
    sigP = varargin{2};
    ei = varargin{3};
end

nullMOI = zeros(NUM_NULL, n);
nullSOI = zeros(NUM_NULL, n);

[inDeg, nullMOI(1, :), nullSOI(1, :)] = scatterInDegVsOutOfIns( ...
    dir_generate_srand(wts), 0);

parfor i = 2:NUM_NULL
    [~, nullMOI(i, :), nullSOI(i, :)] = scatterInDegVsOutOfIns( ...
        dir_generate_srand(wts), 0);
end

finalNullMoiMean = mean(nullMOI);
[finalNullMoiStdL, finalNullMoiStdU] = semistd(nullMOI);
finalNullSoiMean = mean(nullSOI);
[finalNullSoiStdL, finalNullSoiStdU] = semistd(nullSOI);

figure; hold on;
a=7;
b=5;
handles = [];
names = {};
handles(1, 1) = plot(inDeg, finalNullMoiMean, 'ko', 'MarkerSize', b,...
    'LineWidth', 1);
names{1} = 'Null Model Mean';
[~, avgOIM, stdOIM] = scatterInDegVsOutOfIns(wts, 0);
title('Means');
errbar(inDeg, finalNullMoiMean, finalNullMoiStdL, ...
    finalNullMoiStdU, 'k-');
colormap hot;
colormap(flipud(colormap));
c = zeros(1, NUM_NULL);
le = avgOIM < finalNullMoiMean;
ge = avgOIM >= finalNullMoiMean;
c(le) = normcdf(-abs(finalNullMoiMean(le) - avgOIM(le)) ...
    ./ finalNullMoiStdL(le), 0, 1);
c(ge) = normcdf(-abs(avgOIM(ge) - finalNullMoiMean(ge)) ...
    ./ finalNullMoiStdU(ge), 0, 1);
c(c > sigP) = sigP;
scatter(inDeg, avgOIM, a*5, c, 'filled');
hind = 2;
if sum(ei) ~= 0
    handles(1, hind) = plot(inDeg(ei), avgOIM(ei), 'o',  ...
        'MarkerSize', a, 'Color', [.8 .1 .1], 'LineWidth', 1);
    names{hind} = 'Excitatory';
    hind = hind + 1;
end
if sum(~ei) ~= 0
    handles(1, hind) = plot(inDeg(~ei), avgOIM(~ei), 'o',  ...
        'MarkerSize', a, 'Color', [0 .1 .8], 'LineWidth', 1);
    names{hind} = 'Inhibitory';
end
if length(names) == 2
    l = legend(handles, names{1}, names{2});
else
    l = legend(handles, names{1}, names{2}, names{3});
end
set(l, 'Location', 'best');
set(l, 'EdgeColor', [1 1 1]);
h = colorbar;
ylabel(h, 'p-value');
hold off;

figure; hold on;
handles2 = [];
names2 = {};
handles2(1, 1) = plot(inDeg, finalNullSoiMean, 'ko', 'MarkerSize', b ...
    , 'LineWidth', 1);
names2{1} = 'Null Model Mean';
title('Std. Devs');
errbar(inDeg, finalNullSoiMean, finalNullSoiStdL, ...
    finalNullSoiStdU, 'k-');
colormap hot;
colormap(flipud(colormap));
c = zeros(1, NUM_NULL);
le = stdOIM < finalNullSoiMean;
ge = stdOIM >= finalNullSoiMean;
c(le) = normcdf(-abs(finalNullSoiMean(le) - stdOIM(le))  ...
    ./ finalNullSoiStdL(le), 0, 1);
c(ge) = normcdf(-abs(stdOIM(ge) - finalNullSoiMean(ge)) ...
    ./ finalNullSoiStdU(ge), 0 , 1);
c(c > sigP) = sigP;
scatter(inDeg, stdOIM, 5*a, c, 'filled');
hind = 2;
if sum(ei) ~= 0
    handles2(1, hind) = plot(inDeg(ei), stdOIM(ei), 'o',  ...
        'MarkerSize', a, 'Color', [.8 .1 .1], 'LineWidth', 1);
    names2{hind} = 'Excitatory';
    hind = hind + 1;
end
if sum(~ei) ~= 0
    handles2(1, hind) = plot(inDeg(~ei), stdOIM(~ei), 'o',  ...
        'MarkerSize', a, 'Color', [0 .1 .8], 'LineWidth', 1);
    names2{hind} = 'Inhibitory';
end
if length(names2) == 2
    l = legend(handles2, names2{1}, names2{2});
else
    l = legend(handles2, names2{1}, names2{2}, names2{3});
end
set(l, 'Location', 'best');
set(l, 'EdgeColor', [1 1 1]);
h = colorbar;
ylabel(h, 'p-value');
hold off;

end

