
[N, M] = size(wtMat);
ei = (sum(wtMat, 2) > 0)';
stdNumBins = N / 10 + 1;
wts = wtMat(ei, ei);

figure; hist(nonzeros(wts), stdNumBins); title('EPSPs');
figure; hist(PrefFRs, stdNumBins); title('Pref. Firing Rates (spks/s)');
figure; hist(FiringRates, stdNumBins);
title('Firing Rates (spks/s)');

[suIn, mIn, stIn] = statsNZ(wtMat);
[suOut, mOut, stOut] = statsNZ(wtMat');
wtMatRw = dir_generate_srand(wtMat);
[suInRw, mInRw, stInRw] = statsNZ(wtMatRw);
[suOutRw, mOutRw, stOutRw] = statsNZ(wtMatRw');

[kIn, kOut, k] = nodeDegrees(wtMat);
nodeVersatility(wtMat, k);
bins = linspace(0, max(k), stdNumBins);
figure; bar(bins, histc(k, bins), 'histc'); title('Degree');
figure; bar(bins, histc(kIn, bins), 'histc'); title('In Degree');
figure; bar(bins, histc(kOut, bins), 'histc'); title('Out Degree');
allRichClubs(wtMat);
sigScatterInVOut(wtMat, 50);
%neighborHist(wtMat);
%diversity(wtMat, PrefFRs, 1);

N = sum(sum(wtMat(ei, ei) > 0));
p = N / (sum(ei) * (sum(ei) - 1));
biEx = sum(sum((wtMat(ei, ei) > 0) .* (wtMat(ei, ei)' > 0)));
biRatio = biEx / (N * p * p);
wtMatEx = wts;
figure; songetal_fig4;
%dynamicalImportance(wtMat, PrefFRs, 1);

