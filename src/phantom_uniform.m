function [pos, amp] = phantom_uniform(P)
%PHANTOM_UNIFORM  P1: uniform random scatterers in the x-z plane (y = 0).
%   Amplitude is zero-mean Gaussian (sigma = 1) to mimic random impedance.

rng(P.seed, 'twister');
pos = [ P.phX(1) + diff(P.phX) * rand(P.scatN, 1), ...
        zeros(P.scatN, 1), ...
        P.phZ(1) + diff(P.phZ) * rand(P.scatN, 1) ];
amp = randn(P.scatN, 1);
end
