function M = psf_map(P, xg, zg, opt)
%PSF_MAP  复合成像的点扩散函数图：逐点放一个散射子，走完整成像链。
%
%   M = PSF_MAP(P, XG, ZG) 在网格 (XG, ZG) 上逐点放一个单点散射子，跑一遍与正式
%   采集**完全相同**的链路：
%       f2_setup → acquire_pw（发射延时 x1*sin/c）→ das_pw（动态接收孔径）
%       → 11 个角度的射频相加 → 取包络
%   然后统计该点的峰值幅度、相干因子与 -6 dB 宽度。
%
%   这张图回答的问题是：**成像区域内每个点的成像质量是否一致**。
%   幅值差多少 dB、分辨率差多少 mm、复合在整幅上是否都相干，都从这张图读出。
%   它取代了早先那种“把 11 个角度的声场平移后相加”的合成声场图——那幅图的竖条
%   只是“只补了轴上延时”的干涉残差，不是成像里真实存在的加权
%   （见 docs/实验报告_讨论与附录.md 的 B13）。
%
%   输入
%     P    ：params 返回的参数结构体
%     xg   ：横向网格 [m]（行或列都可以），只决定“把散射子摆在哪儿”
%     zg   ：轴向网格 [m]，同上
%     opt  ：可选字段，决定“在每个点周围怎么量宽度”
%              .half  局部窗口的半宽 [m]，默认 1.5e-3（窗口总宽 3 mm）
%              .step  局部窗口的采样间隔 [m]，默认 0.05e-3
%
%   输出（各字段都是 numel(zg) x numel(xg)，单位 SI）
%     PK    ：峰值幅度（绝对值，不做归一化）
%     PdB   ：峰值相对全图最大的 dB 值，20log10(PK/max(PK))
%     CF    ：峰值像素上 11 个角度复值的相干因子 |sum(v)|/sum(|v|)
%     WX    ：峰值处的横向 -6 dB 全宽 [m]
%     WZ    ：峰值处的轴向 -6 dB 全宽 [m]
%     PX,PZ ：峰值在细网格上的实际坐标 [m]，用来检查点位有没有跑偏
%
%   opt.half 与 opt.step 分别管什么
%   -------------------------------
%   每个散射子会临时开一张局部小图：以散射子为中心，x 与 z 都取
%   [x0-half, x0+half] 这一段、按 step 采样（默认 61 x 61 = 3721 个像素）。这张
%   小图就是该点的成像网格，宽度只在它上面量。
%
%   .half 管“窗口装不装得下整个主瓣”：w6 要找峰值两侧降到半幅的两个交点，两个
%   交点都必须落在窗口内，否则返回 NaN。横向主瓣约 0.73 mm 宽（半宽 0.36 mm），
%   1.5 mm 有 4 倍余量；实测把它缩到 0.3 mm，WX 立刻变 NaN，而轴向（瓣宽仅
%   0.195 mm）照常读出 0.196 mm——所以窗口尺寸要按最宽的方向定，不是按最窄的。
%
%   .step 管“量得多细”：它既决定宽度的读数精度，也决定 PX/PZ 的量化步长（PX 只
%   能取 step 的整数倍）。实测（轴上 30 mm 一个点，half 固定 1.5 mm）：
%     step 0.05 mm -> WX 0.730 mm, WZ 0.195 mm, PX 0.000 mm   ← 当前取值
%     step 0.20 mm -> WX 0.866 mm, WZ 0.600 mm, PX -0.100 mm  ← 横向偏大 19 %
%     step 0.40 mm -> WX 0.673 mm, WZ 0.937 mm, PX +0.100 mm  ← 读数失去意义
%   粗采样把窄瓣平滑变宽，读出的已不是真实的 -6 dB 交点。
%
%   代价按点数平方增长：每个点的像素数是 (2*half/step + 1)^2，全图 165 点 x
%   11 次发射合计约 675 万次延时取值，这就是本函数要跑 21–35 s 的原因。所以
%   half 加倍或 step 减半都等于耗时约 4 倍——够用即可，不要随手加密。
%
%   注意：本函数会调用 f2_setup（内部 field_init），重置 Field II 的全局状态；
%   调用前不要持有别的探头句柄，调用后旧句柄也随作废。

if nargin < 4 || isempty(opt), opt = struct(); end
if ~isfield(opt, 'half'), opt.half = 1.5e-3; end
if ~isfield(opt, 'step'), opt.step = 0.05e-3; end

xg = xg(:).';  zg = zg(:);                 % 横向为行、轴向为列，与 das_pw 一致
nx = numel(xg);  nz = numel(zg);
na = P.nAngles;

% 逐角度的时间零点校正（E0 标定）：正式链路里 das_pw 的 dt 就是它，这里必须一致，
% 否则峰值采样点会整体偏移。取不到就退化为 0，不影响幅度（只影响峰值落在哪个采样）。
cf = fullfile(P.dataDir, 'calib.mat');
if exist(cf, 'file'), dt = load(cf, 'dt');  dt = dt.dt(:); else, dt = zeros(na,1); end

[Th, ThR] = f2_setup(P);                   % 与采集完全相同的两个句柄

M = struct('xg', xg, 'zg', zg, 'opt', opt);
M.PK = nan(nz, nx);  M.CF = nan(nz, nx);
M.WX = nan(nz, nx);  M.WZ = nan(nz, nx);
M.PX = nan(nz, nx);  M.PZ = nan(nz, nx);

for iz = 1:nz
    for ix = 1:nx
        x0 = xg(ix);  z0 = zg(iz);
        % 该点的局部成像网格：散射子固定在窗口中心 (x0,z0)，测量只在这张图上做
        xl = x0 + (-opt.half:opt.step:opt.half);     % 默认 61 点，覆盖 ±1.5 mm
        zl = z0 + (-opt.half:opt.step:opt.half);

        % ---- 逐角度：一次发射 + 该角度的 DAS ----------------------------
        A3 = zeros(numel(zl), numel(xl), na);
        for ia = 1:na
            % 散射子只有 1 个，但 acquire_pw 的接口要 N x 3 的坐标与 N x 1 的幅度
            [rf, t0] = acquire_pw(Th, ThR, P, P.angles(ia), [x0 0 z0], 1);
            A3(:,:,ia) = das_pw(rf, t0, P.angles(ia), P, xl, zl, dt(ia));
        end

        % ---- 射频域相加，然后才取包络（顺序不能反，这就是相干复合）--------
        A3 = hilbert(A3);
        Ax = abs(sum(A3, 3));
        [pk, im] = max(Ax(:));
        [i1, i2] = ind2sub(size(Ax), im);

        % ---- 该点上的相干因子与两个方向的宽度 ----------------------------
        v = squeeze(A3(i1, i2, :));
        M.PK(iz,ix) = pk;
        M.CF(iz,ix) = abs(sum(v)) / sum(abs(v));
        M.WX(iz,ix) = w6(xl, Ax(i1,:));
        M.WZ(iz,ix) = w6(zl, Ax(:,i2));
        M.PX(iz,ix) = xl(i2);
        M.PZ(iz,ix) = zl(i1);
    end
    fprintf('    z = %4.0f mm  done\n', zg(iz)*1e3);
end

M.PdB = 20*log10(M.PK / max(M.PK(:)));     % 全图共用同一个参考电平
end

% =========================================================================
function w = w6(g, prof)
%W6  -6 dB 全宽（半幅度阈值），线性插值找两个交点；峰值贴边则返回 NaN。
[pk, i0] = max(prof);
if pk <= 0, w = NaN; return; end
thr = pk/2;
il = i0;  while il > 1 && prof(il) > thr, il = il - 1; end
ir = i0;  while ir < numel(prof) && prof(ir) > thr, ir = ir + 1; end
if il == 1 || ir == numel(prof), w = NaN; return; end     % 单边没有交点：不可测
xl = interp1([prof(il) prof(il+1)], g([il il+1]), thr);
xr = interp1([prof(ir) prof(ir-1)], g([ir ir-1]), thr);
w = abs(xr - xl);
end
