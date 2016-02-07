function [g, gi, go] = degPlot( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Handle inputs, getting or calculating degree values
if length(varargin) == 1
    adjMat = varargin{1} ~= 0;
    if length(size(varargin{1})) > 2
        [n, ~, numSets] = size(adjMat);
    else
        numSets = 1; 
        [n, ~] = size(adjMat);
    end
    k = zeros(numSets, n);
    kIn = zeros(numSets, n);
    kOut = zeros(numSets, n);
    for i = 1:numSets
        [kIn(i,:), kOut(i,:), k(i,:)] = nodeDegrees(adjMat(:,:,i));
    end
elseif length(varargin) == 3
    k = varargin{1};
    kIn = varargin{2};
    kOut = varargin{3};
    [numSets, ~] = size(k);
else
    error('Inappropriate # of input args. 1 input indicates an adjacency matrix (or matrices), 3 indicates kIn, kOut, k.')
end

%% Perform any distribution fitting that requires unaltered k, kIn, kOut
    
pdk = fitdist(nonzeros(k(:)), 'gamma');
pdi = fitdist(nonzeros(kIn(:)), 'gamma');
pdo = fitdist(nonzeros(kOut(:)), 'gamma');
knll = negloglik(pdk);
kinll = negloglik(pdi);
konll = negloglik(pdo);
%pdkn = fitdist(nonzeros(k(:)), 'gp');
%pdin = fitdist(nonzeros(kIn(:)), 'gp');
%pdon = fitdist(nonzeros(kOut(:)), 'gp');
%kdL = negloglik(pdkn);
%kiL = negloglik(pdin);
%koL = negloglik(pdon);


%% Get all histograms, convert to probabilities, and calculate mean and std
[kpdf, kedge] = histcounts(k(:), 'BinMethod', 'fd', 'Normalization', 'pdf'); 
[kipdf, kiedge] = histcounts(kIn(:), 'BinMethod', 'fd', ...
    'Normalization', 'pdf'); 
[kopdf, koedge] = histcounts(kOut(:), 'BinMethod', 'fd', ...
    'Normalization', 'pdf'); 
% Get bin centers...
kbins = diff(kedge)/2+kedge(1:length(kedge)-1);
kibins = diff(kiedge)/2+kiedge(1:length(kiedge)-1);
kobins = diff(koedge)/2+koedge(1:length(koedge)-1);

if numSets > 1
    
    kDeg = zeros(numSets, length(kbins));
    kInDeg = zeros(numSets, length(kibins));
    kOutDeg = zeros(numSets, length(kobins));
    for i = 1:numSets
        kDeg(i, :) = histcounts(k(i, :), kedge, 'Normalization', 'pdf');
        kInDeg(i, :) = histcounts(kIn(i, :), kiedge, ...
            'Normalization', 'pdf');
        kOutDeg(i, :) = histcounts(kOut(i, :), koedge, ...
            'Normalization', 'pdf');
    end

    k_ = mean(kDeg);
    [kL, kU] = semistd(kDeg);
    kIn_ = mean(kInDeg);
    [kInL, kInU] = semistd(kInDeg);
    kOut_ = mean(kOutDeg);
    [kOutL, kOutU] = semistd(kOutDeg);
    

    kL = kL(kbins > 0);
    kU = kU(kbins > 0);
    kInL = kInL(kibins > 0);
    kInU = kInU(kibins > 0);
    kOutL(kOutL < 0) = 0;
    kOutL = kOutL(kobins > 0);
    kOutU = kOutU(kobins > 0);

else
    k_ = kpdf;
    kIn_ = kipdf;
    kOut_ = kopdf;
end
k_ = k_(kbins > 0);
kIn_ = kIn_(kibins > 0);
kOut_ = kOut_(kobins > 0);


[a, b, g] = powerlawfit(kbins(kbins > 0), k_);
[ai, bi, gi] = powerlawfit(kibins(kibins > 0), kIn_);
[ao, bo, go] = powerlawfit(kobins(kobins > 0), kOut_);
[alph, xmin, kdL] = plfit(k_);
[alphi, xmini, kiL] = plfit(kIn_);
[alpho, xmino, koL] = plfit(kOut_);

% Display gofs 
%g.adjrsquare
%gi.adjrsquare
%go.adjrsquare

% Shadder flips out for log y-axis if you don't get rid of zeros...
kNz = find(k_ > 0);
kInNz = find(kIn_ > 0);
kOutNz = find(kOut_ > 0);

%% Plot everything
figure; hold on;
kx = kbins(kbins>0);
kix = kibins(kibins>0);
kox = kobins(kobins>0);

%% Degree
subplot(3, 3, [1 2]); hold on;
if numSets > 1
    posVals = k_ - kL;
    posVals = posVals > 0;
    H = shadedErrorBar(kx(posVals), k_(posVals), [kU(posVals);kL(posVals)], 'k-' , 1);
    set(H.mainLine, 'Color', [0.15 0.15 0.15]);
    set(H.mainLine, 'LineWidth', 2); 
    p1 = H.mainLine;
else
    p1 = plot(kx(kNz), k_(kNz), '-', 'Color', [.3 .3 .3], ...
        'LineWidth', 4);
end
%p2 = plot(kx,pdf(pdkn, kx), '--', 'Color', [0 0 0] , 'LineWidth', 2);
p2 = plot(kx, a*kx.^b, '--', 'Color', [0 0 0] , 'LineWidth', 2);
p3 = plot(kx, pdf(pdk, kx), ':', 'Color', [.5 .5 .5] ,...
    'LineWidth', 3, 'Markersize', 2);
l1 = legend([p1 p2 p3], 'Degree', ['Power-Law (', ...
    num2str(a), ', ', num2str(b), ')'], ['Gamma (', num2str(pdk.a), ...
    ', ', num2str(pdk.b), ')']);
set(l1, 'Interpreter', 'tex');
set(l1, 'Location', 'best');
set(l1, 'EdgeColor', [1 1 1]);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;
subplot(3, 3, 3);
hi = histogram(k, 'BinMethod', 'fd', 'Normalization', 'pdf');
set(hi, 'FaceColor', [0 0 0]);
set(hi, 'EdgeColor', [0 0 0]);
set(hi, 'EdgeAlpha', .6);

%% In Degree
subplot(3, 3, [4 5]); hold on;
if numSets > 1
    posVals = kIn_ - kInL;
    posVals = posVals > 0;
    H = shadedErrorBar(kix(posVals), kIn_(posVals), [kInU(posVals);kInL(posVals)], ...
        'r-', 1); 
    set(H.mainLine, 'LineWidth', 2); 
    p1 = H.mainLine;
else
    p1 = plot(kix(kInNz), kIn_(kInNz), ...
        'r-', 'LineWidth', 4);
end
%p2 = plot(kix, pdf(pdin, kix), '--', 'Color', [.5 0 0] , 'LineWidth', 2);
p2 = plot(kix, ai*kix.^bi, '--', 'Color', [0.5 0 0] , 'LineWidth', 2);
p3 = plot(kix, pdf(pdi, kix), ':', 'Color', [1 .3 .3] , ...
    'LineWidth', 3,'Markersize', 2);
l2 = legend([p1 p2 p3], 'In Degree', ['Power-Law (', ...
    num2str(ai), ', ', num2str(bi)], ['Gamma (' num2str(pdi.a), ...
    ', ', num2str(pdi.b), ')']);
set(l2, 'Interpreter', 'tex');
set(l2, 'Location', 'best');
set(l2, 'EdgeColor', [1 1 1]);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;
subplot(3, 3, 6);
hi = histogram(kIn, 'BinMethod', 'fd', 'Normalization', 'pdf');
set(hi, 'FaceColor', [.5 0 0]);
set(hi, 'EdgeColor', [.5 0 0]);
set(hi, 'EdgeAlpha', .6);

%% Out Degree
subplot(3, 3, [7 8]); hold on;
if numSets > 1
    posVals = kOut_ - kOutL;
    posVals = posVals > 0;
    H = shadedErrorBar(kox(posVals), kOut_(posVals), ...
        [kOutU(posVals);kOutL(posVals)], 'b-', 1);
    set(H.mainLine, 'LineWidth', 2); 
    p1 = H.mainLine;
else
    p1 = plot(kox(kOutNz), kOut_(kOutNz), 'b-',  ...
        'LineWidth', 4);
end
%p2 = plot(kox, pdf(pdon, kox), '--', 'Color', [0 0 .5] , 'LineWidth', 2);
p2 = plot(kox, ao*kox.^bo, '--', 'Color', [0 0 0.5] , 'LineWidth', 2);
p3 = plot(kox, pdf(pdo, kox),':', 'Color', [.3 .3 1], ...
    'LineWidth', 3, 'Markersize', 2);
l3 = legend([p1 p2 p3], 'Out Degree', ['Power-Law (', ...
    num2str(ao), ', ', num2str(bo), ')'], ['Gamma (' num2str(pdo.a), ...
    ', ', num2str(pdo.b), ')']);
set(l3, 'Interpreter', 'tex');
set(l3, 'Location', 'best');
set(l3, 'EdgeColor', [1 1 1]);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;
subplot(3, 3, 9);
hi = histogram(kOut, 'BinMethod', 'fd', 'Normalization', 'pdf');
set(hi, 'FaceColor', [0 0 .5]);
set(hi, 'EdgeColor', [0 0 .5]);
set(hi, 'EdgeAlpha', .6);

end

