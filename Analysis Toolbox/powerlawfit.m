function [ a, b, g ] = powerlawfit(x, y)
    % y = a x^b 
    % Set up fittype and options.
    ft = fittype( 'power1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    % Fit model to data.
    [f, g] = fit( x', y', ft, opts );
    a=f.a;
    b=f.b;
end

