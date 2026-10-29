function run_axisres(root, force)
%RUN_AXISRES  轴上逐点单散射子实验：x = 0 每 2 mm 一个散射子，逐个放入、逐个扫描。
%
%   RUN_AXISRES(ROOT) 在 x = 0、z = 10:2:60 mm 的 26 个深度上，一次只放一个散射子，
%   每个深度都把三种方法（PW / C11 / 常规聚焦）的成像链完整跑一遍，量 -6 dB 横向
%   宽度，存 data/res/axisres.mat。这是图 F10 的数据来源。
%
%   RUN_AXISRES(ROOT, true) 忽略缓存强制重算（约 3.5 min，其中聚焦的 0.5 mm 细扫占 3 min）。
%
%   为什么单独做这个实验
%   ------------------
%   图 F10 原先由点目标体模 P3 轴上那 6 个散射子给出（15、25、30、35、45、55 mm），
%   6 个点画不出曲线形状：焦点附近的拐点、焦外的起伏都落在采样点之间。这里改成
%   2 mm 一点、26 点，曲线才真正连续。数据本身是"逐点重扫"得到的，不是对旧曲线
%   插值，所以也顺带是一次复现性检查。
%
%   它不在 run_core 里
%   -----------------
%   与 run_psfmap 同理：run_core 只管采集与指标，逐点实验由 run_refresh / run_tail
%   调用。冷启动之后必须再跑一次 refresh，否则图 F10 没有数据。
%
%   缓存：结果只依赖参数、深度序列、测量窗口与两个时间零点修正（dt、dt_foc），
%   参数没变就直接用已有文件（见 src/result_stamp.m、src/tag_match.m）。
%
%   输出（data/res/axisres.mat）
%     M.zg/wPW/wC11/wFoc 等，见 src/psf_axis.m；M.tag 是输入指纹（不参与绘图）
%
%   用法：run_axisres

if nargin < 1 || isempty(root)
    root = fileparts(fileparts(mfilename('fullpath')));
end
if nargin < 2, force = []; end
force = force_flag(force);
P = params(root);

% 深度序列：10–60 mm、间隔 2 mm。下界取 10 mm（5–10 mm 两维模型不可靠，见
% params 里 depthMin 的说明），上界取 60 mm（成像深度上限）；30 mm 的焦点正好
% 落在序列上（10 + 2*10 = 30）。
zg  = (10:2:60) * 1e-3;
opt = struct('halfX', 5e-3, 'halfZ', 3e-3, 'amp', 1);   % 与 metrics_resolution 同窗口

% 聚焦方法另做一次 0.5 mm 细扫，只扫聚焦：它的宽度沿深度有约 1.8 mm 周期的起伏
% （硬边发射孔径在焦区的轴向纹波，峰值声压起伏 4–6 dB），2 mm 的 26 个点分辨不出
% 这个尺度。细扫不是插值、也不是"更准的版本"——灰线上的每个点同样是一次独立的
% 单散射子实验，与 2 mm 的点在重合深度上逐位相同；它只是把"两点之间"也测出来，
% 让图上的锯齿看得出是系统的而不是噪声。以它为对照，主结论仍用 2 mm 的 26 个点。
% PW 与 C11 没有焦点，宽度沿深度光滑（实测 2 mm 采样下起伏 <= 0.04 mm），
% 不需要再细分。
zgF = (10:0.5:60) * 1e-3;

out = fullfile(P.dataDir, 'res', 'axisres.mat');
if ~exist(fileparts(out), 'dir'), mkdir(fileparts(out)); end

% 指纹里要含 dt / dt_foc：它们决定峰值落在哪个采样点，宽度会随之微变
cf = fullfile(P.dataDir, 'calib.mat');
if exist(cf, 'file')
    C = load(cf, 'dt'); dt = C.dt(:);
    if isfield(C, 'dt_foc'), dtf = C.dt_foc; else, dtf = 0; end
else
    dt = [];  dtf = [];
end
tag = result_stamp(P, struct('zg', zg, 'zgF', zgF, 'opt', opt, ...
                             'dt', dt, 'dt_foc', dtf));

if ~force && tag_match(out, tag)
    fprintf('\n轴上逐点实验：参数未变，跳过（%s）\n', out);
    return;
end

fprintf('\n===== 轴上逐点单散射子：%d 个深度（%.0f–%.0f mm，间隔 2 mm）=====\n', ...
        numel(zg), zg(1)*1e3, zg(end)*1e3);
t0 = tic;
M = psf_axis(P, zg, opt);

% 第二遍：聚焦方法的 0.5 mm 细扫（只有聚焦；101 个深度，约 3 min——每个深度要
% 发射 58 次，比第一遍的三种方法加起来还费时）
fprintf('\n----- 聚焦方法 0.5 mm 细扫（%d 个深度）-----\n', numel(zgF));
optF = opt;  optF.methods = {'foc'};
Mf = psf_axis(P, zgF, optF);
M.zgFine  = zgF;
M.wFocFine = Mf.wFoc;                       % 细扫只保留宽度，其余字段不需要
M.tag = tag;
save(out, '-struct', 'M', '-v7.3');

% ---- 命令行摘要 --------------------------------------------------------
fprintf('\n深度 (mm)   PW      C11     聚焦\n');
for iz = 1:numel(zg)
    fprintf('  %4.0f     %.3f    %.3f    %.3f\n', ...
            zg(iz)*1e3, M.wPW(iz)*1e3, M.wC11(iz)*1e3, M.wFoc(iz)*1e3);
end
fprintf('\nPW %.3f–%.3f | C11 %.3f–%.3f | 聚焦 %.3f–%.3f mm（2 mm 采样）\n', ...
        min(M.wPW)*1e3, max(M.wPW)*1e3, min(M.wC11)*1e3, max(M.wC11)*1e3, ...
        min(M.wFoc)*1e3, max(M.wFoc)*1e3);
fprintf('聚焦 0.5 mm 细扫 %.3f–%.3f mm（起伏 %.3f mm）\n', ...
        min(M.wFocFine)*1e3, max(M.wFocFine)*1e3, ...
        (max(M.wFocFine)-min(M.wFocFine))*1e3);

% 与表 T-B 的 30 mm 交叉验证：同一套测量口径，新实验应当复现旧读数
i30 = find(abs(zg - 30e-3) < 1e-9, 1);
if ~isempty(i30)
    fprintf('30 mm 交叉验证（T-B: 0.84 / 0.75 / 0.65 mm）: %.3f / %.3f / %.3f mm\n', ...
            M.wPW(i30)*1e3, M.wC11(i30)*1e3, M.wFoc(i30)*1e3);
end
fprintf('saved %s  (%.1f s)\n', out, toc(t0));
end
