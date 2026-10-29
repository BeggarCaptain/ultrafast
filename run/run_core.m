function run_core(mode)
%RUN_CORE  主入口：E0–E5 核心实验。支持断点续跑，并自动识别参数变化。
%
%   在交互式 MATLAB 会话中启动：
%       J = batch(@run_core, 0, {}, 'Profile', 'Processes');
%   或在命令行启动：
%       matlab -batch "addpath('run'); run_core"
%
%   三种模式
%   --------
%   run_core            （= 'auto'）参数与已有数据一致时续跑；不一致直接报错，
%                       并列出改了哪几个参数，绝不静默复用旧数据。
%   run_core('cold')    冷启动：删除全部采集结果（data/pw、data/foc、data/noise、
%                       data/field、data/res、calib.mat、参数指纹）后从头跑，约 60–80 min（视机器负载）。
%   run_core('adopt')   明确接受现有数据：把当前参数登记为数据来源，然后续跑。
%                       只在“确认数据确实由当前参数产生”时使用。
%
%   参数指纹
%   --------
%   采集阶段会写 data/stamp_acq.mat，记录所有影响采集结果的参数（见
%   src/data_stamp.m）与采集算法版本号。改了这些参数（例如 elemPitch、angles、
%   randomSeed、pulseSigmaCycles）后再运行，stamp_guard 会在采集开始前就停下，
%   而不是把旧数据当成新结果。只改显示或分析参数（displayRangeDB、
%   compoundCounts 等）不会触发。
%
%   用 batch 起一个独立 MATLAB worker 有两个好处：交互式 MATLAB 不会被占住；
%   Field II 万一崩了也不会把桌面会话一起带走。
%
%   断点续跑
%   --------
%   每次发射的结果都会在进入下一次之前单独落盘（data/pw/*.mat、data/foc/*.mat），
%   各采集阶段遇到已存在的文件就跳过。任务中断后直接再调用 run_core 即可，它会
%   从第一个缺失的角度/线继续。只想重跑某一次发射，删掉对应文件即可。
%
%   阶段顺序及其理由
%   ----------------
%   1. E0 标定       —— 很快，而且后面所有阶段都依赖它给出的 dt
%   2. P3 点目标     —— 只有 8 个散射子，几分钟就能跑完；先用它验证整条链路
%                       （延时模型、DAS、复合），再去啃两万散射子的体模
%   3. P1 散斑体模   —— 11 个平面波 + 128 条聚焦线（主要耗时）
%   4. REM 囊肿孔洞  —— 从 P1 中挖掉的那部分散射子；P2 = P1 - REM
%   5. 噪声 / 声场   —— 噪声图像、发射声场图
%   6. 分析          —— 复合、指标、表格

if nargin < 1 || isempty(mode)
    mode = 'auto';
end
mode = lower(char(mode));

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);                    % .../ultrafast
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);

% ---- 日志：控制台输出的同时写入 data/logs/core.log --------------------
logdir = fullfile(P.dataDir, 'logs');
if ~exist(logdir, 'dir'), mkdir(logdir); end
diary(fullfile(logdir, 'core.log'));
diary on;

fprintf('=== core run start %s (mode: %s) ===\n', char(datetime('now')), mode);
fprintf('MATLAB %s on %s\n', version, computer);

% ---- 冷启动 / 参数指纹（放在 try 之外：不一致要能让 batch 以非零退出）----
switch mode
    case 'cold'
        pipeline_reset(P, false);
    case 'adopt'
        stamp_guard(P, 'write');
    case 'auto'
        % 交给下面的 stamp_guard 检查
    otherwise
        error('run_core:badMode', '未知模式 "%s"（可选 auto / cold / adopt）', mode);
end
stamp_guard(P);

t0 = tic;

try
    % ---------------------------------------------------------------- E0
    % 只在缺失时标定：它是确定性的，虽然只要几分钟，但里面含 11 次单散射子发射。
    if ~exist(fullfile(P.dataDir,'calib.mat'), 'file')
        stage0_calib(P);
    else
        fprintf('E0: calib.mat exists, skipped\n');
    end

    % ------------------------------------------------- P3（快速验证）
    % 8 个点散射子，很便宜；下面的 QA 行必须打印出约 0.85 mm 与约 0.20 mm
    % （单平面波、−6 dB 口径；0.7 mm 是 11 角度复合后的值，别混用）。
    % 如果数字不对，说明延时模型有问题——此时停下，比三小时后再发现问题好。
    [pp, pa] = phantom_points(P);
    acquire_all_pw(P, 'P3', pp, pa, false);
    acquire_all_focused(P, 'P3', pp, pa, P.psfLineX, false);

    % 在点目标体模上做一次 QA
    Sp = load(fullfile(P.dataDir,'pw',sprintf('P3_a%03d.mat', P.zeroAngleIndex)));
    env = abs(hilbert(double(Sp.img)));     % 提取包络
    m = metrics_resolution(env, P.imageX, P.imageZ, 0, 30e-3);
    fprintf('QA: PW(0 deg) point at 30 mm -> wLat = %.2f mm, wAx = %.2f mm\n', ...
        m.wLat*1e3, m.wAx*1e3);

    % ------------------------------------------------------------- P1 / REM
    % P1 是均匀散斑体模（也用于 SNR）；keepCh=true 是因为 E7 与图 F6
    % 需要原始通道数据。
    [pos, amp] = phantom_uniform(P);
    acquire_all_pw(P, 'P1', pos, amp, true);
    acquire_all_focused(P, 'P1', pos, amp, P.txLineX, true);

    % REM 只含九个囊肿内部的散射子，也就是要从 P1 中减掉、从而得到无回声体模 P2
    % 的那一部分。
    [rp, ra] = phantom_cyst_removed(P, pos, amp);
    fprintf('cyst phantom: %d of %d scatterers removed\n', numel(ra), numel(amp));
    acquire_all_pw(P, 'REM', rp, ra, false);
    acquire_all_focused(P, 'REM', rp, ra, P.txLineX, false);

    % -------------------------------------------------------- 噪声 + 声场
    % 噪声必须在 P1 之后：sigma_e 是用 P1 的信号标定的。
    stage_sigma(P);
    stage_noise(P, 'pw');
    stage_noise(P, 'foc');
    stage_fields(P);

    % ------------------------------------------------------------- 分析
    stage_analysis(P);

    fprintf('=== core run finished in %.1f min ===\n', toc(t0)/60);
catch err
    % 打印可读的调用栈：否则 batch 任务失败时什么信息都没有。
    fprintf(2, 'ERROR: %s\n', err.message);
    for k = 1:numel(err.stack)
        fprintf(2, '  at %s line %d\n', err.stack(k).name, err.stack(k).line);
    end
end

diary off;
end
