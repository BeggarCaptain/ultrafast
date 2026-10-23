function run_refresh()
%RUN_REFRESH  用已有的采集结果重算声场、指标、图与表——完全不重新采集。
%
%   当你只改动了“采集之后”的东西时用它：
%     * 指标定义（例如散斑去趋势的窗口大小）；
%     * 图的布局或色标；
%     * 报告表格的排版；
%     * 复合子集的定义（C11 之类）。
%
%   这里仍然会重跑 stage_fields：合成发射声场只要约 30 s，而且它依赖 calib.mat
%   里的 dt；如果只改了指标，把那次调用删掉即可。
%
%   典型耗时：约一分钟。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
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
    stage_fields(P);          % 发射声场 + 合成声场（调用 calc_hp）
    stage_analysis(P);        % 复合、指标、CSV
    make_figures();           % 图 F1…F13
    make_fig_pointtargets();  % 图 F14（点目标定位验证，只读已存数据）
    make_report_tables();     % 报告用的 Markdown 表格
    fprintf('=== refresh finished in %.1f min ===\n', toc(t0)/60);
catch err
    fprintf(2, 'ERROR: %s\n', err.message);
    for k = 1:numel(err.stack)
        fprintf(2, '  at %s line %d\n', err.stack(k).name, err.stack(k).line);
    end
end
diary off;
end
