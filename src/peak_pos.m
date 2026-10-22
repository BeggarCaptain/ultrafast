function zp = peak_pos(z, env)
%PEAK_POS  Robust axial position of a point-target PSF.
%   The global maximum can jump to a sidelobe; the intensity centroid of the
%   main lobe is unbiased for a symmetric PSF and therefore more stable.

z = z(:);  e = env(:).^2;
[~, i] = max(movmean(env(:), 5));      % smoothed coarse maximum
z0 = z(i);
w  = abs(z - z0) <= 1.2e-3;
zp = sum(z(w).*e(w)) / sum(e(w));
w  = abs(z - zp) <= 1.0e-3;
zp = sum(z(w).*e(w)) / sum(e(w));
end
