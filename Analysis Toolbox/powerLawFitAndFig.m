function [ f, g ] = powerLawFitAndFig( data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [CdfF,CdfX] = ecdf(data,'Function','cdf');
    BinInfo.rule = 1;
    [~,BinEdge] = internal.stats.histbins(data,[],[],BinInfo,CdfF,CdfX);
   
    [y,x] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
    % Set up fittype and options.
    ft = fittype( 'power1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.StartPoint = [0 -1];
    % Fit model to data.
    [f, g] = fit( x', y', ft, opts );
    %pd = fitdist(data, 'gamma');
    f
    g
    figure; hold on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    plot(x, y, 'b-', 'LineWidth', 2);
    %plot(x, pdf(pd, x), 'LineWidth', 2);
    curveBins = logspace(log10(min(x)), log10(max(x)), 20);
    plot(curveBins, f.a.*curveBins.^f.b, ...
        'k--', 'LineWidth', 2);
    hold off;
end

