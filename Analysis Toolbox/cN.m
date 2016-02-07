function [cmmN, cmmNO, cmmNI] = cN (mat)
        [r, c] = size(mat);
        mat = mat ~= 0;
        cmmNO = zeros(r, c);
        cmmNI = zeros(r, c);
        for i = 1:r-1
           cmmNO(i, :) = sum((repmat(mat(i, :), r, 1) .* mat), 2)';
        end
        for i = 1:c-1
           cmmNI(:, i) = sum((repmat(mat(:, i), 1, c) .* mat))';
        end
        cmmNI = triu(cmmNI, 1);
        cmmNO = triu(cmmNO, 1);
        cmmN = cmmNI + cmmNO;
end
