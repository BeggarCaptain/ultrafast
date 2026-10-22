function [rf, t0, d] = acquire_focused(Th, P, xf, pos, amp)
%ACQUIRE_FOCUSED  One focused transmit, returns the 128-channel raw RF.
%   The transmit aperture is a rectangle of width P.txApodW centred on XF
%   (52 elements, clipped at the array edges), i.e. F# = 1.75 on transmit.
%   The firing delays are set explicitly (not with xdc_focus) so that the very
%   same vector can be removed again in the DAS: Field II adds the firing delay
%   to each channel's time axis (measured in E0).

F = [xf 0 P.txFocus];
d = (norm(F) - sqrt((P.xe - xf).^2 + P.txFocus^2)) / P.c;

xdc_center_focus(Th, [0 0 0]);
xdc_focus_times(Th, 1, d);

apo = double(abs(P.xe - xf) <= P.txApodW/2 + 1e-12);
xdc_apodization(Th, 1, apo);

[rf, t0] = calc_scat_multi(Th, Th, pos, amp);
end
