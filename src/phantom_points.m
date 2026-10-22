function [pos, amp] = phantom_points(P)
%PHANTOM_POINTS  P3: eight strong point scatterers.

pos = [P.ptX(:), zeros(numel(P.ptX), 1), P.ptZ(:)];
amp = P.ptAmp * ones(numel(P.ptX), 1);
end
