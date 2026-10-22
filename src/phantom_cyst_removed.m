function [pos, amp, sel] = phantom_cyst_removed(P, pos0, amp0)
%PHANTOM_CYST_REMOVED  Scatterers that fall inside the nine anechoic cysts.
%   The nine cysts are the Cartesian product of P.cystX and P.cystZ (3 x 3).
%   Because Field II is linear in scatterer amplitude, the P2 channel data is
%   exactly  P1 - REM  (see design section 4.4).

C = cyst_centres(P);
sel = false(size(pos0, 1), 1);
for k = 1:size(C, 1)
    d2 = (pos0(:,1) - C(k,1)).^2 + (pos0(:,3) - C(k,2)).^2;
    sel = sel | (d2 <= P.cystR^2);
end
pos = pos0(sel, :);
amp = amp0(sel);
end
