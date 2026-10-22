function s = pulse_stats(rf, fs)
%PULSE_STATS  Centre frequency and -6 dB bandwidth of the two-way pulse.
%   Uses the strongest channel of a single-point acquisition.

[~, ich] = max(max(abs(rf), [], 1));
w   = rf(:, ich) .* tukeywin(size(rf, 1), 0.1);
nfft = 8 * numel(w);
S   = abs(fft(w, nfft));
f   = (0:nfft-1)' * fs / nfft;

half = f < fs/2;
f = f(half);  S = S(half);
S = S / max(S);

idx = find(S >= 10^(-3/20));
s.f1 = f(idx(1));
s.f2 = f(idx(end));
s.fc = (s.f1 + s.f2) / 2;
s.bw = (s.f2 - s.f1) / s.fc;
s.f  = f;
s.S  = 20*log10(max(S, 1e-12));
s.wav = w;
s.ich = ich;
end
