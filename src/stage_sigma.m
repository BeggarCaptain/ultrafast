function stage_sigma(P)
%STAGE_SIGMA  Calibrate the channel noise level and validate the analytic
%   noise-propagation formula by Monte-Carlo (design 4.5).
%   target: single plane wave, z = 30 mm -> SNR = P.snrPW30 dB.

fprintf('\n===== noise calibration =====\n');
C = load(fullfile(P.dataDir,'calib.mat'));
S = load(fullfile(P.dataDir,'pw','P1_a036.mat'));   % alpha = 0
img = double(S.img);  w2 = double(S.w2);

[XX, ZZ] = meshgrid(P.x, P.z);
roi = abs(XX) <= 5e-3 & ZZ >= 25e-3 & ZZ <= 35e-3;

Ps = mean(img(roi).^2);
W  = mean(w2(roi));
sigma_e = sqrt(Ps / (10^(P.snrPW30/10) * W));
fprintf('signal power %.3e, mean sum(w^2) %.1f -> sigma_e %.3e\n', Ps, W, sigma_e);

% ---------------------------------------------------- Monte-Carlo validation
xr = (-5:0.1711:5)*1e-3;  zr = (25:0.1:35)*1e-3;
Nt = size(S.ch, 1);
[~, w2r] = das_pw(zeros(Nt,P.N,'single'), S.t0, 0, P, xr, zr, S.dt);
acc = zeros(numel(zr), numel(xr));  acc2 = acc;
rng(11,'twister');
for r = 1:P.nMC
    nz = sigma_e * randn(Nt, P.N);
    nr = das_pw(nz, S.t0, 0, P, xr, zr, S.dt);
    acc  = acc  + nr;
    acc2 = acc2 + nr.^2;
end
emp = sqrt(max(acc2/P.nMC - (acc/P.nMC).^2, 0));
ana = sigma_e * sqrt(w2r);
ratio = emp ./ ana;
fprintf('Monte-Carlo (%d reps): empirical/analytic noise ratio = %.3f +- %.3f\n', ...
    P.nMC, mean(ratio(:)), std(ratio(:)));

C.sigma_e = sigma_e;
C.mcRatio = mean(ratio(:));
C.mcRatioStd = std(ratio(:));
C.Ps30 = Ps;  C.W30 = W;
save(fullfile(P.dataDir,'calib.mat'), '-struct', 'C');
fprintf('sigma_e saved to calib.mat\n');
end
