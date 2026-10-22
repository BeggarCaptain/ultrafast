function run_refresh()
%RUN_REFRESH  Recompute field maps, metrics, figures and tables from the
%   existing acquisitions (no Field II re-acquisition).

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'run'));
P = params(root);

logdir = fullfile(P.dataDir, 'logs');
if ~exist(logdir, 'dir'), mkdir(logdir); end
diary(fullfile(logdir, 'refresh.log')); diary on;
fprintf('=== refresh start %s ===\n', char(datetime('now')));
t0 = tic;
try
    stage_fields(P);
    stage_analysis(P);
    make_figures();
    make_report_tables();
    fprintf('=== refresh finished in %.1f min ===\n', toc(t0)/60);
catch err
    fprintf(2, 'ERROR: %s\n', err.message);
    for k = 1:numel(err.stack)
        fprintf(2, '  at %s line %d\n', err.stack(k).name, err.stack(k).line);
    end
end
diary off;
end
