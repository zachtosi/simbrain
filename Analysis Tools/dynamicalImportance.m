function [ dynImp ] = dynamicalImportance( wtMat, vars, show )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
        adj = single(wtMat ~= 0);
        maxEig = max(real(eig(adj)));
        [N, ~] = size(wtMat);
        maxEig = maxEig(1);
        
        dynImp = zeros(1, N);
        
        parfor i = 1:N
            fracEig = max(real(eig(adj(1:N ~= i, 1:N ~= i))));
            dynImp(i) = fracEig / maxEig;
        end
        dynImp = 1 - dynImp;

        if show == 1
            
           figure; hist(dynImp, N/10 + 1); title('Dynamical Importance');
           figure; scatter(vars, dynImp); 
           percentiles = zeros(1, N);
           for i = 1:N
              percentiles(i) = 100 * sum(vars < i) / N; 
           end
           figure; scatter(dynImp, percentiles); title('Percentiles'); 
        end

end

