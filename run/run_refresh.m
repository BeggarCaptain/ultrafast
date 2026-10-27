function run_refresh()
%RUN_REFRESH  用已有的采集结果重算声场、指标、图与表——完全不重新采集。
%
%   当你只改动了“采集之后”的东西时用它：
%     * 指标定义（例如分辨率宽度取 −10 dB 还是 −6 dB）；
%     * 图的布局或色标；
%     * 报告表格的排版；
%     * 复合子集的定义（C11 之类）。
%
%   这里仍然会重跑 stage_fields：11 个角度的发射声场只要约 4 s，而且它依赖
%   Field II 的阵列配置；如果只改了指标，把那次调用删掉即可。
%
%   也会重跑 run_psfmap：图 F5 要在成像区域内逐点成像，约 35 s。它跑的是真实
%   的 Field II 链路（每点 11 次发射），是这条流程里最慢的一步；只改图的样式时
%   同样可以去掉它，data/res/psfmap.mat 里已有结果。
%
%   典型耗时：约两分钟。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'run'));
P = params(root);
stamp_guard(P);        % 只重算分析也必须参数一致，否则等于拿旧数据出报告

logdir = fullfile(P.dataDir, 'logs');
if ~exist(logdir, 'dir'), mkdir(logdir); end
diary(fullfile(logdir, 'refresh.log')); diary on;
fprintf('=== refresh start %s ===\n', char(datetime('now')));
t0 = tic;
try
    stage_fields(P);          % 各角度的发射声场（调用 calc_hp）
    stage_analysis(P);        % 复合、指标、CSV
    run_psfmap(root);         % 成像区域逐点成像（图 F5 的数据）
    make_figures();           % 图 F1…F12（含 F5）
    make_fig_coverage(root);  % 图 F14（成像加权与覆盖范围）
    make_fig_pointtargets();  % 图 F13（点目标定位验证，只读已存数据）
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
