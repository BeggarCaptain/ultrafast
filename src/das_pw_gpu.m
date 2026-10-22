function img = das_pw_gpu(rf, t0, alpha, P, dt)
%DAS_PW_GPU  GPU version of DAS_PW (E7 only; no noise-weight output).
%   RF must be a gpuArray.  Identical delay model to das_pw.m.

z = P.z(:);  x = P.x(:).';
[Nt, Nch] = size(rf);
Nz = numel(z);  Nx = numel(x);

img = gpuArray.zeros(Nz, Nx, 'single');
ca = cos(alpha);  sa = sin(alpha);
off = gpuArray((0:Nch-1) * Nt);
dxe = gpuArray(P.xe);

for k = 1:Nx
    xc  = x(k);
    dxc = dxe - xc;
    hw  = z / (2*P.Fnum);
    m   = abs(dxc) <= hw;

    tau = (z*ca + xc*sa)/P.c + sqrt(z.^2 + dxc.^2)/P.c + dxc*sa/P.c + dt;
    pos = (tau - t0) * P.fs + 1;

    i0 = floor(pos);
    fr = pos - i0;
    m  = m & (i0 >= 1) & (i0 < Nt);
    i0 = min(max(i0, 1), Nt-1);

    v0 = rf(i0 + off);
    v1 = rf(i0 + 1 + off);
    img(:, k) = sum(((1-fr).*m).*v0 + (fr.*m).*v1, 2);
end
end
