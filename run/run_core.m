function run_core()
%RUN_CORE  Master entry point: E0-E5 (core experiments).  Resume-safe.
%   Launch from MATLAB with
%       J = batch(@run_core, 0, {}, 'Profile', 'Processes');
%   The function is resume-safe: every angle/line is stored separately, so a
%   restarted run continues where the previous one stopped.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);

logdir = fullfile(P.dataDir, 'logs');
if ~exist(logdir, 'dir'), mkdir(logdir); end
diary(fullfile(logdir, 'core.log'));
diary on;

fprintf('=== core run start %s ===\n', char(datetime('now')));
fprintf('MATLAB %s on %s\n', version, computer);
t0 = tic;

try
    % ---------------------------------------------------------------- E0
    if ~exist(fullfile(P.dataDir,'calib.mat'), 'file')
        stage0_calib(P);
    else
        fprintf('E0: calib.mat exists, skipped\n');
    end

    % ------------------------------------------------- P3 (fast, sanity)
    [pp, pa] = phantom_points(P);
    acquire_all_pw(P, 'P3', pp, pa, false);
    acquire_all_focused(P, 'P3', pp, pa, P.xpsf, false);

    % quick QA on the point-target phantom
    Sp = load(fullfile(P.dataDir,'pw','P3_a036.mat'));
    env = abs(hilbert(double(Sp.img)));
    m = metrics_resolution(env, P.x, P.z, 0, 30e-3);
    fprintf('QA: PW(0 deg) point at 30 mm -> w10 = %.2f mm, w6 = %.2f mm\n', ...
        m.w10*1e3, m.w6*1e3);

    % ------------------------------------------------------------- P1 / REM
    [pos, amp] = phantom_uniform(P);
    acquire_all_pw(P, 'P1', pos, amp, true);
    acquire_all_focused(P, 'P1', pos, amp, P.xf, true);

    [rp, ra] = phantom_cyst_removed(P, pos, amp);
    fprintf('cyst phantom: %d of %d scatterers removed\n', numel(ra), numel(amp));
    acquire_all_pw(P, 'REM', rp, ra, false);
    acquire_all_focused(P, 'REM', rp, ra, P.xf, false);

    % -------------------------------------------------------- noise + fields
    stage_sigma(P);
    stage_noise(P, 'pw');
    stage_noise(P, 'foc');
    stage_fields(P);

    % ------------------------------------------------------------- analysis
    stage_analysis(P);

    fprintf('=== core run finished in %.1f min ===\n', toc(t0)/60);
catch err
    fprintf(2, 'ERROR: %s\n', err.message);
    for k = 1:numel(err.stack)
        fprintf(2, '  at %s line %d\n', err.stack(k).name, err.stack(k).line);
    end
end

diary off;
end
