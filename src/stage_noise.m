function stage_noise(P, kind)
%STAGE_NOISE  White channel noise -> beamformed noise images.
%   KIND = 'pw'  : data/noise/N1_a###.mat, one per plane-wave angle
%   KIND = 'foc' : data/noise/NF_l###.mat, one per focused line
%   The same realisations are used for P1 and P2 (P2 = P1 - REM), as in a real
%   acquisition.

C = load(fullfile(P.dataDir, 'calib.mat'));
outdir = fullfile(P.dataDir, 'noise');
if ~exist(outdir, 'dir'), mkdir(outdir); end
sig = C.sigma_e;
Th  = f2_setup(P);
rng(7, 'twister');

if strcmp(kind, 'pw')
    for i = 1:P.nAngles
        fs_ = fullfile(outdir, sprintf('N1_a%03d.mat', i));
        if exist(fs_, 'file'), continue; end
        S0 = load(fullfile(P.dataDir, 'pw', sprintf('P1_a%03d.mat', i)));
        Nt = size(S0.img, 1);
        if isfield(S0, 'ch'), Nt = size(S0.ch, 1); end
        nz = sig * randn(Nt, P.N);
        [img, ~] = das_pw(nz, S0.t0, S0.alpha, P, P.x, P.z, S0.dt);
        S = struct('img', single(img), 'tag', 'N1', 'alpha', S0.alpha);
        save(fs_, '-struct', 'S');
        if mod(i,10)==0, fprintf('[noise-pw] %d/%d\n', i, P.nAngles); end
    end
else
    for k = 1:P.txLines
        fs_ = fullfile(outdir, sprintf('NF_l%03d.mat', k));
        if exist(fs_, 'file'), continue; end
        S0 = load(fullfile(P.dataDir, 'foc', sprintf('P1_l%03d.mat', k)));
        % channel-support length: only the centre line stores 'ch', so fall
        % back to the plane-wave record (same phantom, same array, same range)
        if isfield(S0, 'Nt')
            Nt = S0.Nt;
        elseif isfield(S0, 'ch')
            Nt = size(S0.ch, 1);
        else
            Sp = load(fullfile(P.dataDir, 'pw', 'P1_a001.mat'));
            Nt = size(Sp.ch, 1);
        end
        nz = sig * randn(Nt, P.N);
        [img, ~] = das_focused(nz, S0.t0, S0.xf, P, P.z, S0.dt);
        S = struct('img', single(img), 'tag', 'NF', 'xf', S0.xf);
        save(fs_, '-struct', 'S');
        if mod(k,20)==0, fprintf('[noise-foc] %d/%d\n', k, P.txLines); end
    end
end
end
