function run_tail(mode)
%RUN_TAIL  在 P1 阶段之后收尾：囊肿体模、噪声、声场、分析、绘图与表格。
%   支持断点续跑，并检查参数指纹。
%
%   RUN_TAIL          默认 'auto'：参数与已有数据不一致时报错
%   RUN_TAIL('adopt') 接受现有数据，把当前参数登记为数据来源（有风险）
%
%   冷启动请用 run_core('cold')——它会把 P1 的数据也一起清掉。
%
%   为什么单独做一个入口
%   --------------------
%   P1 采集（11 个平面波 + 128 条聚焦线，两万散射子）是最贵的部分。这部分一旦
%   落盘，其余都可以便宜地重跑：囊肿体模（REM）只需要两千多个散射子，分析加
%   画图也就一分钟。单独拆出 run_tail 后，改一个指标定义或改一张图的布局，都
%   不必重跑数小时的采集。
%
%   它同时也是 run_core 在 P1 之后中断时的恢复路径：REM 采集自身断点续跑，
%   噪声/声场/分析三个阶段要么可续跑、要么是确定性的。
%
%   三个阶段为什么按这个顺序
%   ----------------------
%   REM → 噪声 → 声场 → 分析。REM 必须用自己的散射子集合重新跑一遍发射；噪声
%   依赖 run_core 里已经标定好的 sigma_e（P1 之后才算得出来），所以不能提到前面；
%   声场与分析都只读磁盘上的文件，放在最后。每一步都是"文件存在就跳过"，所以
%   中途失败再跑一次不会重做已经完成的部分。
%
%   为什么 REM 要用与 P1 相同的随机种子重建
%   --------------------------------------
%   phantom_uniform 内部会重设种子，因此这里拿到的两万个散射子与 P1 采集时完全
%   相同；囊肿孔洞（REM）是从这一批里挑出来的那两千多个。P2 = P1 − REM 之所以
%   成立，前提就是两者用的是同一批散射子——否则相减会把两个不同体模的散斑搅在
%   一起，出来的不是"无回声囊肿"而是噪声。
%
%   它不画 F13：点目标定位图由 make_fig_pointtargets 单独生成，走这条路时要补
%   一次调用（README 的复现小节写明了）。

if nargin < 1 || isempty(mode)
    mode = 'auto';
end
mode = lower(char(mode));

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
fprintf('=== tail run start %s (mode: %s) ===\n', char(datetime('now')), mode);

switch mode
    case 'adopt'
        stamp_guard(P, 'write');
    case 'auto'
        stamp_guard(P);
    otherwise
        error('run_tail:badMode', ...
            'run_tail 只支持 auto / adopt；冷启动请用 run_core(''cold'')');
end
preflight_tail(P);      % 缺 P1 采集结果或 sigma_e 时给出可读的提示
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
    acquire_all_focused(P, 'REM', rp, ra, P.txLineX, false);

    % ---- 噪声、声场、指标、图、表 --------------------------------------
    stage_noise(P, 'pw');
    stage_noise(P, 'foc');
    stage_fields(P);
    stage_analysis(P);
    run_psfmap(root);         % 成像区域逐点成像（图 F5 的数据，21–35 s；命中指纹则跳过）
    run_axisres(root);        % 轴上逐点单散射子实验（图 F10 的数据；命中指纹则跳过）
    make_figures();           % 图 F1…F12（含 F5）
    make_fig_coverage(root);  % 图 F14（成像加权与覆盖范围）
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

% =========================================================================
function preflight_tail(P)
%PREFLIGHT_TAIL  检查 run_tail 的前置条件：P1 的采集结果 + 已标定的 sigma_e。
%
%   run_tail 只负责“P1 之后”的部分：REM 采集、噪声、声场、指标、图、表。它假定
%   P1（散斑体模的 11 个平面波 + 128 条聚焦线）以及 calib.mat 里的 sigma_e 已经
%   存在。缺这些文件时直接跑会在 stage_noise 里报一个很难懂的 load 错误，所以
%   这里先查一遍，给出“先跑 run_core”的提示。

need = cell(1, P.nAngles + numel(P.txLineX) + 1);
for i = 1:P.nAngles
    need{i} = fullfile(P.dataDir, 'pw', sprintf('P1_a%03d.mat', i));
end
for k = 1:numel(P.txLineX)
    need{P.nAngles + k} = fullfile(P.dataDir, 'foc', sprintf('P1_l%03d.mat', k));
end
need{end} = fullfile(P.dataDir, 'calib.mat');

miss = need(~cellfun(@(f) exist(f, 'file') == 2, need));
if ~isempty(miss)
    error('run_tail:noP1', ['run_tail 只做 P1 之后的收尾，但缺少 %d 个前置文件，' ...
        '例如：\n  %s\n请先跑 run_core（它会从缺失的那个角度/线继续）。'], ...
        numel(miss), miss{1});
end

C = load(fullfile(P.dataDir, 'calib.mat'));
if ~isfield(C, 'sigma_e')
    error('run_tail:noSigma', ['calib.mat 里没有 sigma_e —— 它由 stage_sigma 在 P1 ' ...
        '之后标定。\n请先跑一次 run_core（P1 已存在时会直接跳到噪声标定，约 1 分钟）。']);
end
end
