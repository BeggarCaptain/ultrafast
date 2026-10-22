function [rf, t0] = acquire_pw(Th, P, alpha, pos, amp)
%ACQUIRE_PW  One plane-wave transmit, returns the 128-channel raw RF.
%   rf : [Nt x N] double, t0 : time of the first sample [s].

d = P.xe * sin(alpha) / P.c;         % element firing delays (row, may be < 0)
xdc_center_focus(Th, [0 0 0]);
xdc_focus_times(Th, 0, d);
xdc_apodization(Th, 1, ones(1, P.N));   % full aperture, no apodisation

[rf, t0] = calc_scat_multi(Th, Th, pos, amp);
end
