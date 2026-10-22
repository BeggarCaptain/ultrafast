function stage0_calib(P)
%STAGE0_CALIB  E0: pulse, delay model, time origin, apertures, interpolation.
%   Writes data/calib.mat with everything the later stages need.

fprintf('\n===== E0 calibration =====\n');
Th = f2_setup(P);

% ------------------------------------------------------------------ (1) pulse
ptpos = [0 0 30e-3];  ptamp = 1;
[rf1, t01] = acquire_pw(Th, P, 0, ptpos, ptamp);
fprintf('single-scatterer record: %d samples x %d channels, t0 = %.3f us\n', ...
    size(rf1,1), size(rf1,2), t01*1e6);

ps = pulse_stats(rf1, P.fs);
fprintf('pulse: fc = %.3f MHz, -6 dB BW = %.1f %%\n', ps.fc/1e6, 100*ps.bw);

% ------------------------------------------------------------- (2) tilt check
xchk = (-15:1:15)*1e-3;
pchk = [xchk(:), zeros(numel(xchk),1), 40e-3*ones(numel(xchk),1)];
tilt = zeros(1,3);  slope = zeros(1,3);
alphas_chk = [0 10 16.5]*pi/180;
for q = 1:3
    d = P.xe*sin(alphas_chk(q))/P.c;
    xdc_center_focus(Th,[0 0 0]); xdc_focus_times(Th, 0, d);
    xdc_apodization(Th, 1, ones(1,P.N));
    [hp, th0] = calc_hp(Th, pchk);
    tt = th0 + (0:size(hp,1)-1)'/P.fs;
    [~, im] = max(abs(hp), [], 1);
    tp = tt(im);
    cf = polyfit(xchk(:), tp(:), 1);
    slope(q) = cf(1);
    if abs(sin(alphas_chk(q))) > 1e-9
        tilt(q) = cf(1) / (sin(alphas_chk(q))/P.c);
    else
        tilt(q) = NaN;
    end
    fprintf('alpha=%4.1f deg: measured slope %.4f us/mm, theory %.4f us/mm, ratio %.3f\n', ...
        alphas_chk(q)*180/pi, cf(1)*1e6, sin(alphas_chk(q))/P.c*1e6, tilt(q));
end

% --------------------------------------- (3) time-origin calibration per angle
% With dt = 0 the point at 30 mm appears at z_obs; the model delay must be
% increased by delta = 2 (z_obs - 30 mm)/c so that the image peak lands on the
% true position (see design 4.2, time origin).
zcal = (26:0.05:34)'*1e-3;
dt = zeros(1, P.nAngles);
for i = 1:P.nAngles
    [rf, t0] = acquire_pw(Th, P, P.angles(i), ptpos, ptamp);
    img = das_pw(rf, t0, P.angles(i), P, 0, zcal, 0);
    dt(i) = 2*(peak_pos(zcal, abs(hilbert(img))) - 30e-3)/P.c;
end
% reject isolated lobe-picking outliers (should be a smooth function of alpha)
for i = 1:P.nAngles
    lo = max(1, i-2);  hi = min(P.nAngles, i+2);
    med = median(dt(lo:hi));
    if abs(dt(i) - med) > 0.5e-6
        fprintf('  dt outlier at angle %d: %.3f -> %.3f us\n', i, dt(i)*1e6, med*1e6);
        dt(i) = med;
    end
end
fprintf('dt(alpha): min %.3f us, max %.3f us, mean %.3f us, spread %.3f us\n', ...
    min(dt)*1e6, max(dt)*1e6, mean(dt)*1e6, (max(dt)-min(dt))*1e6);

% focused reference
[rfF, t0F] = acquire_focused(Th, P, 0, ptpos, ptamp);
imgF = das_focused(rfF, t0F, 0, P, zcal, 0);
dt_foc = 2*(peak_pos(zcal, abs(hilbert(imgF))) - 30e-3)/P.c;
fprintf('dt(focused, xf=0) = %.3f us\n', dt_foc*1e6);

% verification: the calibrated peak must land on 30.00 mm for every angle
zver = zeros(1,4);  aver = [1 19 36 71];
for q = 1:4
    i = aver(q);
    [rf, t0] = acquire_pw(Th, P, P.angles(i), ptpos, ptamp);
    img = das_pw(rf, t0, P.angles(i), P, 0, zcal, dt(i));
    zver(q) = peak_pos(zcal, abs(hilbert(img)));
end
fprintf('calibrated peak depth (angles 1/19/36/71): %s mm\n', ...
    mat2str(round(zver*1e3, 3)));

% --------------------------------------------- (4) point-target PSF after calib
xg = (-5:0.05:5)*1e-3;  zg = (26:0.05:34)'*1e-3;   % column vector!
[rf, t0] = acquire_pw(Th, P, 0, ptpos, ptamp);
img = das_pw(rf, t0, 0, P, xg, zg, dt(P.iAng==0));
env = to_env(img);
[~, imax] = max(env(:));  [iz, ix] = ind2sub(size(env), imax);
lat = env(iz, :)';  axp = env(:, ix);
fprintf('PW alpha=0 point: peak at (%.2f, %.2f) mm\n', xg(ix)*1e3, zg(iz)*1e3);
fprintf('  -10 dB lateral width %.2f mm, -6 dB axial width %.2f mm\n', ...
    crossing_width(xg, lat, 10^(-10/20))*1e3, crossing_width(zg, axp, 0.5)*1e3);

% lateral PSF of the focused line requires several lines: done later in E4.

% --------------------------------------------------------- (5) aperture counts
for zc = [10 30 60]*1e-3
    fprintf('aperture at z=%2.0f mm: %3d elements (%.1f mm)\n', zc*1e3, ...
        sum(abs(P.xe) <= zc/(2*P.Fnum)), 2*zc/(2*P.Fnum)*1e3);
end

% ------------------------------------------------------- (6) linear vs spline
[rf, t0] = acquire_pw(Th, P, 0, ptpos, ptamp);
zc = P.z;
imgLin = das_pw(rf, t0, 0, P, 0, zc, dt(P.iAng==0));
imgSp  = zeros(numel(zc),1);
tt = t0 + (0:size(rf,1)-1)'/P.fs;
for c = 1:size(rf,2)
    tau = zc/P.c + sqrt(zc.^2 + P.xe(c)^2)/P.c + dt(P.iAng==0);
    m = abs(P.xe(c)) <= zc/(2*P.Fnum);
    v = interp1(tt, rf(:,c), tau, 'spline', 0);
    imgSp = imgSp + m.*v;
end
d1 = abs(hilbert(imgLin)); d2 = abs(hilbert(imgSp));
band = zc > 20e-3 & zc < 50e-3;
err = 20*log10(max(abs(d1(band)-d2(band))) / max(d1(band)));
fprintf('linear vs spline beamforming: max difference %.1f dB (20-50 mm)\n', err);

% ------------------------------------------------------------ (7) timing test
[pos, amp] = phantom_uniform(P);
tic; [rf20, t020] = acquire_pw(Th, P, 0, pos, amp); tAcq = toc;
fprintf('calc_scat_multi, %d scatterers, %d ch: %.1f s -> record %d samples\n', ...
    numel(amp), P.N, tAcq, size(rf20,1));

% ---------------------------------------------------------------- save result
S = struct();
S.matlabVersion = version;
S.computer      = computer;
S.fieldIIPath   = which('field_init');
S.pulse         = rmfield(ps, {'f','S','wav'});
S.pulseF        = ps.f;  S.pulseS = ps.S;
S.pulseWav      = ps.wav;
S.t01           = t01;
S.tiltRatio     = tilt;
S.tiltSlope     = slope;
S.alphaChk      = alphas_chk;
S.dt            = dt;
S.dt_foc        = dt_foc;
S.pwPeak        = [xg(ix) zg(iz)];
S.pwWidth10     = crossing_width(xg, lat, 10^(-10/20));
S.pwWidth6      = crossing_width(zg, axp, 0.5);
S.interpErrDB   = err;
S.tAcq20k       = tAcq;
S.recordSamples = size(rf20,1);
S.zc            = zc;
S.imgLin        = imgLin;
S.imgSpline     = imgSp;

if ~exist(P.dataDir,'dir'), mkdir(P.dataDir); end
save(fullfile(P.dataDir,'calib.mat'), '-struct', 'S');
fprintf('saved %s\n', fullfile(P.dataDir,'calib.mat'));
end

function e = to_env(x)
e = abs(hilbert(x));
end
