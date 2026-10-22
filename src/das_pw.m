function [img, w2] = das_pw(rf, t0, alpha, P, x, z, dt)
%DAS_PW  Delay-and-sum beamforming of one plane-wave transmit.
%   IMG  : [numel(z) x numel(x)] real RF after receive dynamic focusing.
%   W2   : sum of squared interpolation weights, used for noise propagation.
%   DT   : global time-origin correction [s], calibrated in E0.
%
%   Delay model (design eq. 1-4):
%       tau = (z cos(alpha) + x sin(alpha))/c + sqrt(z^2+(x-xe)^2)/c + dt
%   Receive aperture: rectangular window, |xe-x| <= z/(2F).

if nargin < 7 || isempty(dt), dt = 0; end

z = z(:);  x = x(:).';          % force column depth / row lateral
[Nt, Nch] = size(rf);
Nz = numel(z);  Nx = numel(x);
img = zeros(Nz, Nx);
w2  = zeros(Nz, Nx);

ca = cos(alpha);  sa = sin(alpha);
off = (0:Nch-1) * Nt;
dxe = P.xe;                              % 1 x N
ax  = abs(dxe);

for k = 1:Nx
    xc = x(k);
    axc = abs(dxe - xc);                 % 1 x N
    hw  = z / (2*P.Fnum);                % Nz x 1
    m   = axc <= hw;                     % Nz x N  (implicit expansion)

    % Field II returns each channel with its own transmit firing delay added to
    % the time axis (verified in E0: residual 0.01 us with the term, 1.2 us
    % without).  The delay must therefore be part of the imaging model.
    tau = (z*ca + xc*sa)/P.c + sqrt(z.^2 + (dxe - xc).^2)/P.c ...
          + dxe*sa/P.c + dt;
    pos = (tau - t0) * P.fs + 1;         % 1-based sample index

    i0 = floor(pos);
    fr = pos - i0;
    m  = m & (i0 >= 1) & (i0 < Nt);

    i0 = min(max(i0, 1), Nt-1);
    v0 = rf(i0 + off);
    v1 = rf(i0 + 1 + off);

    w0 = (1 - fr) .* m;
    w1 = fr .* m;

    img(:, k) = sum(w0.*v0 + w1.*v1, 2);
    w2(:, k)  = sum(w0.^2 + w1.^2, 2);
end
end
