function run_tail()
%RUN_TAIL  Finish the core run after the P1 stage: cyst phantom, noise,
%   field maps, analysis, figures and tables.  Resume-safe.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'run'));
P = params(root);

logdir = fullfile(P.dataDir, 'logs');
if ~exist(logdir, 'dir'), mkdir(logdir); end
diary(fullfile(logdir, 'tail.log'));
diary on;
fprintf('=== tail run start %s ===\n', char(datetime('now')));
t0 = tic;

try
    [pos, amp] = phantom_uniform(P);
    [rp, ra] = phantom_cyst_removed(P, pos, amp);
    fprintf('cyst phantom: %d of %d scatterers removed (9 cysts)\n', numel(ra), numel(amp));

    acquire_all_pw(P, 'REM', rp, ra, false);
    acquire_all_focused(P, 'REM', rp, ra, P.xf, false);

    stage_noise(P, 'pw');
    stage_noise(P, 'foc');
    stage_fields(P);
    stage_analysis(P);
    make_figures();
    make_report_tables();

    fprintf('=== tail run finished in %.1f min ===\n', toc(t0)/60);
catch err
    fprintf(2, 'ERROR: %s\n', err.message);
    for k = 1:numel(err.stack)
        fprintf(2, '  at %s line %d\n', err.stack(k).name, err.stack(k).line);
    end
end

diary off;
end
