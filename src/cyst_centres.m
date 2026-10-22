function C = cyst_centres(P)
%CYST_CENTRES  All 9 cyst centres: the Cartesian product of P.cystX and P.cystZ.

[XX, ZZ] = meshgrid(P.cystX, P.cystZ);
C = [XX(:), ZZ(:)];
end
