function [ nn, pn ] = commonNeigh( adjMat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    RW = 100;
    [m, n] = size(adjMat);
    if m~=n
        error('Input matrix is not square');
    end
    [neigh, neighO, neighI] = cN(adjMat);
    ne= zeros(1, m);
    neO = zeros(1, m);
    neI = zeros(1, m);
    am = triu((adjMat+adjMat')~=0, 1);
    pN = 0:(m-1);
    for k = 1:m         
        ne(1, k) = sum(sum((neigh==k-1) .* am))/ ...
            sum(sum((neigh == k-1)));
        neO(1, k) = sum(sum((neighO==k-1) .* am))/ ...
            sum(sum((neighO == k-1)));
        neI(1, k) = sum(sum((neighI==k-1) .* am))/ ...
            sum(sum((neighI == k-1)));
    end
    
    
    neighRw = zeros(RW, m);
    neighORw = zeros(RW, m);
    neighIRw = zeros(RW, m);
    
    parfor j = 1:RW
        rwMat = dir_generate_srand(adjMat);
        [rw, rwO, rwI] = cN(rwMat);
        rwMat = triu((rwMat + rwMat')~=0, 1);
        for k = 1:m         
            neighRw(j, k) = sum(sum((rw==k-1) .* rwMat))/ ...
                sum(sum((rw == k-1)));
            neighORw(j, k) = sum(sum((rwO==k-1) .* rwMat))/ ...
                sum(sum((rwO == k-1)));
            neighIRw(j, k) = sum(sum((rwI==k-1) .* rwMat))/ ...
                sum(sum((rwI == k-1)));
            if isnan(neighRw(j, k))
               neighRw(j, k) = 0; 
            end
            if isnan(neighORw(j, k))
               neighORw(j, k) = 0; 
            end
            if isnan(neighIRw(j, k))
               neighIRw(j, k) = 0; 
            end
        end
        
    end
  
    for k = 1:3
        figure; hold on;
        if k==1
            nrw = neighRw;
            n = ne;
            col = 'k';
        elseif k==2
            nrw = neighORw;
            n = neO;
            col = 'b';
        else
            nrw = neighIRw;
            n = neI;
            col = 'r';
        end
        [stnL, stnU] = semistd(nrw);
        mn = mean(nrw);
        nnans = ~isnan(stnL .* stnU .*mn);
        nnans = nnans & n > 0;
        errbar(pN(nnans), mn(nnans), stnL(nnans), stnU(nnans), col);
        scatter(pN(nnans), mn(nnans), col);
        plot(pN(nnans), n(nnans), col);
        hold off;
    end

end

