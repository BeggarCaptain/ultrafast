function idx = subset_indices(P)
%SUBSET_INDICES  Angle-index sets for the 7 compounding levels (design 4.3).
%   Returned as 1-based indices into P.angles.

I = cell(1, numel(P.subsetN));
for k = 1:numel(P.subsetN)
    switch P.subsetN(k)
        case 1,  ii = 0;
        case 5,  ii = [0 8 -8 16 -16];
        case 9,  ii = [0 4 -4 8 -8 12 -12 16 -16];
        case 17, ii = [0 2:2:16, -2:-2:-16];
        case 33, ii = [0 1:16, -1:-16];
        case 45, ii = [0 1:22, -1:-22];
        case 71, ii = -35:35;
        otherwise, error('unknown subset %d', P.subsetN(k));
    end
    [~, loc] = ismember(ii, P.iAng);
    I{k} = sort(loc);
end
idx = I;
end
