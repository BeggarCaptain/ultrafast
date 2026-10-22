function m = metrics_resolution(env, x, z, xp, zp)
%METRICS_RESOLUTION  -10 dB lateral and -6 dB axial width of one point target.
%   The profiles are up-sampled 8x (shape-preserving) before the crossings are
%   located, otherwise the 0.05 mm axial grid quantises a 0.2 mm wide envelope.

ix = find(abs(x - xp) <= 5e-3);
iz = find(abs(z - zp) <= 3e-3);
E  = double(env(iz, ix));
[~, im] = max(E(:));
[r, c] = ind2sub(size(E), im);

xu = linspace(x(ix(1)), x(ix(end)), 8*numel(ix));
zu = linspace(z(iz(1)), z(iz(end)), 8*numel(iz));
lat = interp1(x(ix).', E(r, :).', xu, 'pchip');
axp = interp1(z(iz),    E(:, c),    zu, 'pchip');

[~, iL] = max(lat);  m.xpk = xu(iL);
[~, iA] = max(axp);  m.zpk = zu(iA);
m.w10 = crossing_width(xu, lat, 10^(-10/20));
m.w6  = crossing_width(zu, axp, 0.5);
m.peak = max(lat);
end
