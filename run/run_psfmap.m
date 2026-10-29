function run_psfmap(root, force)
%RUN_PSFMAP  在成像区域内逐点放散射子，跑完整成像链，得到复合成像的 PSF 图。
%
%   RUN_PSFMAP(ROOT) 在 x = -10:2:10 mm、z = 15:2.5:50 mm 的网格上逐点放一个
%   单散射子，每个点都走一遍正式链路（11 个角度发射 → 128 通道 → DAS → 射频域
%   相加 → 包络），统计峰值幅度、相干因子 CF、横向/轴向 -6 dB 宽度，存
%   data/res/psfmap.mat。
%
%   RUN_PSFMAP(ROOT, true) 忽略缓存强制重算（约 21–35 s，视机器负载）。
%
%   这张图不在 run_core 里
%   ---------------------
%   run_core 只做“采集 + 指标”，它产出 data/pw、data/foc、data/noise 与
%   data/field/maps.mat；逐点成像是一次独立的验证实验（165 个点 × 11 次发射），
%   由 run_refresh / run_tail 调用。所以冷启动之后必须再跑一次 run_refresh，
%   否则图 F5 没有数据（README“复现”一节）。
%
%   缓存：结果只依赖参数与网格，参数没变就直接用已有 psfmap.mat，跳过这 35 s。
%   指纹见 src/result_stamp.m，比较见 src/tag_match.m。
%
%   网格为什么取这一块
%   ------------------
%   11 束平面波要“都照到”才算完整复合。孔径半宽 21.12 mm、最大偏转 10°，在深度
%   z 处 11 束的公共区域是 |x| <= 21.12 - z*tan(10°) mm：z = 50 mm 时还有
%   ±12.3 mm，z = 60 mm 时剩 ±10.5 mm。所以取 |x| <= 10 mm、z = 15..50 mm 这一块，
%   整块都落在“11 束全照到”的范围内，与图 F14 里画的虚线框一致。
%
%   输出（data/res/psfmap.mat）
%     M.PK/PdB/CF/WX/WZ/PX/PZ，均为 numel(zg) x numel(xg)，见 src/psf_map.m
%     M.tag：算这份结果用的输入指纹（不参与绘图）
%
%   用法：run_psfmap

if nargin < 1 || isempty(root)
    root = fileparts(fileparts(mfilename('fullpath')));
end
if nargin < 2, force = []; end
force = force_flag(force);
P = params(root);

xg = (-10:2:10) * 1e-3;          % 横向 11 点
zg = (15:2.5:50) * 1e-3;         % 轴向 15 点
opt = struct('half', 1.5e-3, 'step', 0.05e-3);   % 量宽度的细网格，见 psf_map

out = fullfile(P.dataDir, 'res', 'psfmap.mat');
if ~exist(fileparts(out), 'dir'), mkdir(fileparts(out)); end

% 指纹要含三样不在 P 里的东西：网格、细网格定义、以及 E0 标定的逐角度时间零点
% dt——dt 变了峰值落在哪个采样会变，宽度随之微变。
cf = fullfile(P.dataDir, 'calib.mat');
if exist(cf, 'file'), C = load(cf, 'dt'); dt = C.dt(:); else, dt = []; end
tag = result_stamp(P, struct('xg', xg, 'zg', zg, 'half', opt.half, ...
                             'step', opt.step, 'dt', dt));

if ~force && tag_match(out, tag)
    fprintf('\nPSF map: 参数未变，跳过（%s）\n', out);
    return;
end

fprintf('\n===== compounded PSF map: %d x %d 点，%d 个角度 =====\n', ...
        numel(zg), numel(xg), P.nAngles);
tic;
M = psf_map(P, xg, zg, opt);
M.tag = tag;                     % 与结果一起落盘，下次据此判断能否跳过
save(out, '-struct', 'M', '-v7.3');

% ---- 命令行摘要：成像区域内“各个点”到底有多均匀 -------------------------
fprintf('\n峰值幅度 %.1f .. %.1f dB（跨度 %.2f dB，全图归一化）\n', ...
        min(M.PdB(:)), max(M.PdB(:)), max(M.PdB(:)) - min(M.PdB(:)));
for iz = 1:numel(zg)
    fprintf('  z = %4.1f mm : 横向 %+.1f dB | 中值 %+6.2f dB | CF %.4f | WX %.2f mm\n', ...
            zg(iz)*1e3, max(M.PdB(iz,:)) - min(M.PdB(iz,:)), median(M.PdB(iz,:)), ...
            median(M.CF(iz,:)), median(M.WX(iz,:))*1e3);
end
fprintf('CF %.4f .. %.4f | WX %.2f .. %.2f mm | WZ %.2f .. %.2f mm\n', ...
        min(M.CF(:)), max(M.CF(:)), min(M.WX(:))*1e3, max(M.WX(:))*1e3, ...
        min(M.WZ(:))*1e3, max(M.WZ(:))*1e3);
fprintf('saved %s  (%.1f s)\n', out, toc);
end
