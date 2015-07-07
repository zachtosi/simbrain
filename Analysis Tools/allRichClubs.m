function allRichClubs ( wtMat )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    figure;
    subplot(3, 2, 1)
    richClub(wtMat ~= 0);
    title('In-Rich Club: All Wts')
    
    subplot(3, 2, 2)
    richClub(wtMat' ~= 0);
    title('Out-Rich Club: All Wts')
    
    subplot(3, 2, 3)
    richClub(wtMat > 0);
    title('In-Rich-Club: Ex Only')
    
    subplot(3, 2, 4)
    richClub(wtMat' > 0);
    title('Out-Rich-Club: Ex Only')
    
    wtsSort = sort(nonzeros(abs(wtMat)), 'descend');
    cutOff = wtsSort(int32(numel(wtsSort) / 10));
    
    subplot(3, 2, 5)
    richClub(abs(wtMat) > cutOff);
    title('In-Rich Club: Top 10%')
    
    subplot(3, 2, 6)
    richClub(abs(wtMat)' > cutOff);
    title('Out-Rich Club: Top 10%')
    
end

