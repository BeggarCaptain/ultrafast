function [img, w2] = das_focused(rf, t0, xf, P, z, dt)
%DAS_FOCUSED  Delay-and-sum for one focused B-mode line.
%   Delay model (design eq. 7): tau = (z + sqrt(z^2+(xf-xe)^2))/c + dt
%   Returns IMG [numel(z) x 1] and the weight-square sum W2.

if nargin < 6 || isempty(dt), dt = 0; end

z = z(:);                       % force column depth
[Nt, Nch] = size(rf);
Nz = numel(z);
dxe = abs(P.xe - xf);
hw  = z / (2*P.Fnum);
m   = dxe <= hw;                            % Nz x N
off = (0:Nch-1) * Nt;

tau = (z + sqrt(z.^2 + dxe.^2))/P.c + dt ...
      + (sqrt(xf^2 + P.txFocus^2) - sqrt((xf - P.xe).^2 + P.txFocus^2))/P.c;
pos = (tau - t0) * P.fs + 1;
i0  = floor(pos);
fr  = pos - i0;
m   = m & (i0 >= 1) & (i0 < Nt);

i0 = min(max(i0, 1), Nt-1);
v0 = rf(i0 + off);
v1 = rf(i0 + 1 + off);
w0 = (1 - fr) .* m;
w1 = fr .* m;

img = sum(w0.*v0 + w1.*v1, 2);
w2  = sum(w0.^2 + w1.^2, 2);
end
