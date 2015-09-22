
NUM_NULL = 100;
nullMV = zeros(NUM_NULL, numel(PrefFRs));
for i = 1:NUM_NULL
   nullMV(i, :) = nodeVersatility(dir_generate_srand(wtMat), k, 0); 
end

fig = figure; hold;
scatter(k, mean(nullMV), 'bo');
nullStV = std(nullMV);
nullMV = mean(nullMV);
for i = 1:numel(PrefFRs)
   plot([k(i) k(i)], [nullMV(i)+nullStV(i) nullMV(i)-nullStV(i)], 'b');
end
vers =  nodeVersatility(wtMat, k, 0);
%plotyy(k, vers, k, abs(vers - nullMV) ./ nullStV, @scatter, @stem);
a = 25;
c = abs(vers - nullMV) ./ nullStV;
scatter(k, vers, a, c, 'filled');
hold;