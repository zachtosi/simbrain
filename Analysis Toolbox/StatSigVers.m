function [ vers, stdVer, nullMV, nullStVL, nullStVU ] = StatSigVers( ...
    wtMat, versVar, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    narginchk(2, 5);
    NUM_NULL = 100;
    sigP = .1;
    plt = 1;
    if nargin >= 3
       plt = varargin{1}; 
    end
    if nargin >= 4
        NUM_NULL = varargin{2};
    end
    if nargin == 5
        sigP = varargin{3};
    end
        
    ei = (sum(wtMat, 2)>0)';
    n = numel(versVar);
    nullMV = zeros(NUM_NULL, n);
    parfor i = 1:NUM_NULL
       nullMV(i, :) = nodeVersatility(dir_generate_srand(wtMat), versVar, 0); 
    end
    [nullStVL, nullStVU] = semistd(nullMV);
    nullMV = mean(nullMV);
    vers =  nodeVersatility(wtMat, versVar, 0);
    le = vers < nullMV;
    ge = vers >= nullMV;
    stdVer = zeros(1, n);
    stdVer(le) = (nullMV(le) - vers(le)) ./ nullStVL(le);
    stdVer(ge) = (vers(ge) - nullMV(ge)) ./ nullStVU(ge);

    if plt == 1
        figure;
        subplot(1, 2, 1);
        hold on;
        colormap parula;
        colormap(flipud(colormap));
        a = 50;
        b = 15;
        scatter(versVar, nullMV, b, 'ko');
        errbar(versVar, nullMV, nullStVL, nullStVU, 'k-');  
        c = normcdf(-stdVer, 0, 1);
        c(c>sigP) = sigP;
        scatter(versVar, vers, a, c, 'filled');
        h = colorbar;
        ylabel(h, 'p-value');
        %scatter(versVar(ei), vers(ei), a, [.8 .1 .1]);
        %scatter(versVar(~ei), vers(~ei), a, [0 .1 .8]);
        hold off;
        subplot(1, 2, 2);
        [~, edges] = histcounts(vers, 'Normalization', 'pdf',  ...
            'BinMethod', 'fd');
        hold on;
        hv = histogram(vers, edges, 'Normalization', 'pdf');
        set(hv, 'FaceAlpha', .5);
        set(hv, 'FaceColor', [.1 .4 .6]);
        set(hv, 'EdgeAlpha', 0);
        hv = histogram(nullMV, edges, 'Normalization', 'pdf');
        set(hv, 'FaceAlpha', .5);
        set(hv, 'FaceColor', [0 0 0]);
        set(hv, 'EdgeAlpha', 0);
        
    end

end

