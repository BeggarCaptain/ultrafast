function run_refresh(force)
%RUN_REFRESH  用已有的采集结果重算声场、指标、图与表——**不重跑任何采集**。
%
%   当你只改动了“采集之后”的东西时用它：
%     * 指标定义（例如分辨率宽度取 −10 dB 还是 −6 dB）；
%     * 图的布局或色标；
%     * 报告表格的排版；
%     * 复合子集的定义（C11 之类）。
%
%   “不重跑采集”指什么
%   ------------------
%   data/pw、data/foc、data/noise 一个文件都不动：不调用任何 acquire_*，也不重
%   标定 sigma_e。开始处只做一次 stamp_guard：参数与现有采集数据不一致就直接报错
%   退出，因为那种情况下旧数据本来就不该再用（此时应跑 run_core 或 run_core('cold')）。
%
%   但它**会**重算三样 Field II 的东西，它们是计算而不是采集
%   ------------------------------------------------------
%   * stage_fields：11 个角度的发射声压图，约 4 s，图 F3/F4/F14 要用；
%   * run_psfmap  ：成像区域内 165 个点各跑 11 次发射的点扩散函数图，21–35 s，
%                   图 F5 要用。**run_core 不做这一步**（它只管采集与指标），所以
%                   冷启动之后必须有一次 refresh 才能出 F5；
%   * run_axisres ：轴上 26 个深度的逐点单散射子实验（另加聚焦的 0.5 mm 细扫），
%                   首次约 3.5 min，图 F10 要用。同样只属于 refresh，不属于 run_core。
%
%   这三步都存了一份输入指纹：参数、网格与计算逻辑没变就直接用磁盘上的结果，
%   合计 25–40 s 会省掉（逐点实验首次跑则要几分钟）。想强制重算（换了 Field II
%   版本、要核对结果）时：
%       run_refresh('force')
%
%   典型耗时：全部命中缓存 10 s；强制重算 5–7 min（其中逐点实验占大头）。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
if nargin < 1, force = []; end
force = force_flag(force);
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
    stage_fields(P, force);   % 各角度的发射声场（calc_hp；有缓存）
    stage_analysis(P);        % 复合、指标、CSV
    run_psfmap(root, force);  % 成像区域逐点成像（图 F5 的数据；有缓存）
    run_axisres(root, force); % 轴上逐点单散射子实验（图 F10 的数据；有缓存）
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
