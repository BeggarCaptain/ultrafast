function s = metrics_speckle(env, x, z, xc, zc, half)
%METRICS_SPECKLE  First-order speckle statistics and correlation lengths.
%   The intensity is divided by a locally smoothed version first: without this
%   the 10 mm ROI straddles a several-dB brightness trend (no TGC in the
%   simulation) and mu/sigma drops to ~0.5 even for fully developed speckle.

if nargin < 6 || isempty(half), half = 5e-3; end

I = double(env).^2;

% ---- remove the slow depth/lateral brightness trend ----------------------
kz = max(3, 2*round(1e-3/(z(2)-z(1))) + 1);      % ~2 mm axial box
kx = max(3, 2*round(1e-3/(x(2)-x(1))) + 1);      % ~2 mm lateral box
B  = imfilter(I, ones(kz,kx)/(kz*kx), 'replicate');
In = I ./ max(B, realmin);

[XX, ZZ] = meshgrid(x, z);
sel = abs(XX - xc) <= half & abs(ZZ - zc) <= half;
v = In(sel);

s.mu    = mean(v);
s.sigma = std(v);
s.ratio = s.mu / s.sigma;          % intensity: exponential speckle -> 1.00
s.ratioA = mean(sqrt(v)) / std(sqrt(v));   % envelope: Rayleigh -> 1.91
s.n     = numel(v);

% ---- correlation lengths from the mean profile along each axis -----------
rows = any(sel, 2);  cols = any(sel, 1);
profZ = mean(In(rows, cols), 2);
profX = mean(In(rows, cols), 1);

s.corrZ = corr_len(z(rows), profZ);
s.corrX = corr_len(x(cols), profX);

% ---- one-dimensional intensity histogram (normalised) ---------------------
s.histEdges = linspace(0, 4, 41);
s.hist = histcounts(v / s.mu, s.histEdges);
s.hist = s.hist / sum(s.hist);
end

function w = corr_len(v, p)
v = v(:);  p = p(:) - mean(p);
ac = xcorr(p, 'coeff');
ac = ac(numel(p):end);            % zero lag .. end
vv = (0:numel(ac)-1)' * (v(2)-v(1));
if ac(1) < 0.5
    w = NaN; return;
end
i = find(ac < 0.5, 1);
if isempty(i) || i == 1
    w = NaN; return;
end
w = interp1([ac(i-1); ac(i)], [vv(i-1); vv(i)], 0.5, 'linear', 'extrap');
end
