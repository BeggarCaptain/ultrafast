function w = crossing_width(v, p, lvl)
%CROSSING_WIDTH  Width of a 1-D profile at LVL (linear, relative to peak).
%   V : coordinate vector, P : profile, LVL : e.g. 10^(-10/20) for -10 dB.

v = v(:);  p = p(:);
[~, ip] = max(p);
thr = lvl * p(ip);
n = numel(p);

iL = ip;
while iL > 1 && p(iL) > thr, iL = iL - 1; end
if iL == ip, w = NaN; return; end
vL = interp1([p(iL); p(iL+1)], [v(iL); v(iL+1)], thr, 'linear', 'extrap');

iR = ip;
while iR < n && p(iR) > thr, iR = iR + 1; end
if iR == ip, w = NaN; return; end
vR = interp1([p(iR-1); p(iR)], [v(iR-1); v(iR)], thr, 'linear', 'extrap');

w = vR - vL;
end
