function run_tail()
%RUN_TAIL  在 P1 阶段之后收尾：囊肿体模、噪声、声场、分析、绘图与表格。
%   支持断点续跑。
%
%   为什么单独做一个入口
%   --------------------
%   P1 采集（71 个平面波 + 128 条聚焦线，两万散射子）是最贵的部分。这部分一旦
%   落盘，其余都可以便宜地重跑：囊肿体模（REM）只需要两千多个散射子，分析加
%   画图也就一分钟。单独拆出 run_tail 后，改一个指标定义或改一张图的布局，都
%   不必重跑数小时的采集。
%
%   它同时也是 run_core 在 P1 之后中断时的恢复路径：REM 采集自身断点续跑，
%   噪声/声场/分析三个阶段要么可续跑、要么是确定性的。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
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
    % ---- 必须用与 P1 相同的种子重建囊肿体模 ----------------------------
    % phantom_uniform 内部会重设随机数种子，因此得到的散射子集合与 P1 采集时
    % 完全一致；这样“P1 减去 REM”才是严格成立的相减。
    [pos, amp] = phantom_uniform(P);
    [rp, ra] = phantom_cyst_removed(P, pos, amp);
    fprintf('cyst phantom: %d of %d scatterers removed (9 cysts)\n', numel(ra), numel(amp));

    % ---- 只对被移除的散射子做采集 --------------------------------------
    acquire_all_pw(P, 'REM', rp, ra, false);
    acquire_all_focused(P, 'REM', rp, ra, P.xf, false);

    % ---- 噪声、声场、指标、图、表 --------------------------------------
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
