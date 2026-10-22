function run_dense_snr()
%RUN_DENSE_SNR  Supplementary check: does a denser scatterer phantom recover
%   the ideal 10log10(n) SNR gain of coherent compounding?
%   Same array, same noise level, 20 000 scatterers confined to 20 x 30 mm
%   (33 /mm^2, ~4x the density of P1).  Only plane-wave transmits are needed.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'run'));
P = params(root);
C = load(fullfile(P.dataDir, 'calib.mat'));

outdir = fullfile(P.dataDir, 'dense');
if ~exist(outdir, 'dir'), mkdir(outdir); end
logdir = fullfile(P.dataDir, 'logs');
diary(fullfile(logdir, 'dense.log')); diary on;
fprintf('=== dense phantom SNR check %s ===\n', char(datetime('now')));
t0 = tic;

try
    rng(P.seed + 1, 'twister');
    n = P.scatN;
    pos = [ -10e-3 + 20e-3*rand(n,1), zeros(n,1), 15e-3 + 30e-3*rand(n,1) ];
    amp = randn(n,1);
    fprintf('density %.1f /mm2\n', n/(20*30));

    Th = f2_setup(P);
    for i = 1:P.nAngles
        f = fullfile(outdir, sprintf('D_a%03d.mat', i));
        if exist(f, 'file'), continue; end
        [rf, t00] = acquire_pw(Th, P, P.angles(i), pos, amp);
        [img, w2] = das_pw(rf, t00, P.angles(i), P, P.x, P.z, C.dt(i));
        nz = C.sigma_e * randn(size(rf));
        niv = das_pw(nz, t00, P.angles(i), P, P.x, P.z, C.dt(i));
        S = struct('img', single(img), 'nz', single(niv), 'w2', single(w2), ...
                   't0', t00, 'dt', C.dt(i));
        save(f, '-struct', 'S');
        if mod(i,10) == 0
            fprintf('  angle %2d/%d  elapsed %.1f min\n', i, P.nAngles, toc(t0)/60);
        end
    end

    % ---- SNR gain per subset -------------------------------------------
    sets = subset_indices(P);
    zdep = (15:5:45)*1e-3;
    gain = zeros(1, numel(P.subsetN));
    [zz, s0] = snr_curve_(load(fullfile(outdir,'D_a036.mat')), P);
    ref = interp1(zz, s0, 30e-3);
    for q = 1:numel(P.subsetN)
        S = zeros(numel(P.z), numel(P.x), 'single');
        N = zeros(numel(P.z), numel(P.x), 'single');
        for i = sets{q}
            T = load(fullfile(outdir, sprintf('D_a%03d.mat', i)));
            S = S + T.img;  N = N + T.nz;
        end
        [zz, s] = snr_curve_rm(S, N, P);
        gain(q) = interp1(zz, s, 30e-3) - ref;
        fprintf('  n=%2d: relative SNR at 30 mm = %5.2f dB (theory %.2f)\n', ...
            P.subsetN(q), gain(q), 10*log10(P.subsetN(q)));
    end
    save(fullfile(P.dataDir,'res','dense_snr.mat'), 'P', 'gain', 'zdep');
    fprintf('=== dense check done in %.1f min ===\n', toc(t0)/60);
catch err
    fprintf(2, 'ERROR: %s\n', err.message);
end
diary off;
end

function [zz, snr] = snr_curve_(S, P)
[zz, snr] = snr_curve_rm(S.img, S.nz, P);
end

function [zz, snr] = snr_curve_rm(rf, nrf, P)
ps = double(rf).^2;  pn = double(nrf).^2;
h1 = ones(round(5e-3/P.dz),1); h1 = h1/sum(h1);
h2 = ones(1, round(5e-3/P.dx));  h2 = h2/sum(h2);
ps = imfilter(imfilter(ps, h1, 'replicate'), h2, 'replicate');
pn = imfilter(imfilter(pn, h1, 'replicate'), h2, 'replicate');
S  = 10*log10(ps ./ max(pn, realmin));
zz = P.z;
snr = median(S, 2);
end
