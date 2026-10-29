function M = psf_axis(P, zg, opt)
%PSF_AXIS  轴上逐深度单散射子实验：一次只放一个散射子，每放一次扫一遍。
%
%   M = PSF_AXIS(P, ZG) 在 x = 0、深度 ZG(k) 放**一个**散射子，把整条成像链重跑
%   一遍，再在峰值所在行上量 -6 dB 横向宽度。三种方法各得一条曲线：
%     * PW  ：0° 单次平面波发射（1 次发射）
%     * C11 ：11 个角度发射，射频域相干相加（11 次发射）
%     * 聚焦：逐线发射、逐线重建（本函数只扫中心 ±5 mm 的 58 条线，见下）
%
%   它解决什么问题：图 F10 原来只有 6 个深度点
%   ------------------------------------------
%   旧 F10 的数据来自点目标体模 P3 轴上那 6 个散射子（15、25、30、35、45、55 mm）。
%   6 个点画不出曲线的形状：焦点附近的拐点、焦外的起伏都落在采样点之间，读图的人
%   只能看到"折线"。这里改成 x = 0 上每 2 mm 一个散射子、共 26 个（10–60 mm），
%   得到一条真正连续的曲线；代价是每个深度都要重新发射、重新成像一遍。
%
%   为什么一次只放一个散射子
%   ----------------------
%   两方面的原因。其一，多目标体模里相邻目标的旁瓣会互相叠加——PW 在浅层的
%   "翅膀"尤其明显（图 F13 上排左幅），翅膀交叠处测到的已经不是单个点的 PSF。
%   其二，逐点放、逐点扫，每个深度都是一次独立的实验，测到的宽度只由该深度的
%   声场与孔径决定，不掺入邻居的贡献。
%
%   为什么横向用 P.imageX（0.171 mm 的显示网格）
%   -------------------------------------------
%   为了让这张图和表 T-B 用同一把尺子。T-B 与旧 F10 都是在 0.171 mm 的显示网格上、
%   用 metrics_resolution 的 ±5 mm 窗口与 8 倍 pchip 上采样量的 -6 dB 全宽；这里
%   照搬同一套口径，30 mm 处读数为 PW 0.846 / C11 0.751 / 聚焦 0.645 mm，与 T-B
%   的 0.84 / 0.75 / 0.65 mm 一致——这就是新实验对旧实验的交叉验证。
%   注意 F5 的 WX 读数是 0.71–0.74 mm，那是同一块主瓣在 0.05 mm 局部细网格上的值：
%   采样更细、pchip 插值的平滑量更小，读数略窄。两套数字不矛盾，各自内部自洽。
%
%   为什么聚焦只扫 ±5 mm 的 58 条线
%   -----------------------------
%   聚焦 B 模式每条线的数据只由**它自己那一次发射**决定：Field II 里各次发射的
%   散射计算互相独立，没有跨发射的耦合。所以扫满 256 条线，在中心这 58 条线上
%   得到的结果与这里逐位相同，多花的只是 13 倍时间。宽度只由主瓣决定（0.65–0.90 mm），
%   ±5 mm 的窗口对主瓣有 5 倍以上的余量。
%
%   输入
%     P    ：params 返回的参数结构体
%     zg   ：深度向量 [m]（列或行都可以），只决定把散射子摆在哪
%     opt  ：可选字段
%              .halfX 横向窗口半宽 [m]，默认 5e-3（与 metrics_resolution 一致）
%              .halfZ 轴向窗口半宽 [m]，默认 3e-3（同上）
%              .amp   散射子幅度，默认 1（只量宽度，幅度取多少不影响读数）
%              .methods 要跑哪几种方法，默认三种全跑；只想扫聚焦曲线时传
%                       {'foc'}，省掉 PW 与 C11 的发射时间（每点约 0.9 s）
%
%   输出（各字段都是 numel(zg) x 1）
%     wPW,wC11,wFoc ：三种方法的 -6 dB 横向全宽 [m]（NaN 表示交点落在窗口外）
%     xPW,xC11,xFoc ：实测峰值横向位置 [m]（检查有没有跑偏）
%     pPW,pC11,pFoc ：峰值包络（绝对刻度，只作合理性检查）
%     xg, zs 规则、opt、tag 之外的其余字段见下；xg 是本函数实际用的线坐标
%
%   时间口径
%   --------
%   平面波用 calib.dt(ia)（每个角度一个时间零点），聚焦用 calib.dt_foc；这两者就是
%   正式链路里 das_pw / das_focused 用的值。少加它们，峰值会整体挪到相邻采样点，
%   宽度随之微变（0.02 mm 量级）。
%
%   代价：每个深度 0.6–1.1 s（PW 约 0.4 s、C11 约 0.5 s、聚焦约 0.25 s，随机器负载浮动），
%   26 个深度合计 16–30 s；聚焦的 0.5 mm 细扫每个深度要发 58 次、约 2 s，101 个深度约 3 min。
%
%   注意：本函数调用 f2_setup 会重置 Field II 的全局状态；调用前不要持有别的探头
%   句柄，调用后旧句柄也随作废。

if nargin < 3 || isempty(opt), opt = struct(); end
if ~isfield(opt, 'halfX'), opt.halfX = 5e-3; end   % 与 metrics_resolution 的窗口一致
if ~isfield(opt, 'halfZ'), opt.halfZ = 3e-3; end
if ~isfield(opt, 'amp'),   opt.amp   = 1;    end   % 宽度与幅度无关
if ~isfield(opt, 'methods'), opt.methods = {'pw','c11','foc'}; end
doPW  = any(strcmpi(opt.methods, 'pw'));
doC11 = any(strcmpi(opt.methods, 'c11'));
doFoc = any(strcmpi(opt.methods, 'foc'));

zg = zg(:).';                                     % 横向为行、轴向为列
nz = numel(zg);

% 横向线坐标：显示网格里 |x| <= 5 mm 的那些线（0.171 mm 间距，共 58 条）。
% 直接用 P.imageX 的子集，而不是另造一张网格——这样读数与 T-B 同源。
xg = P.imageX(abs(P.imageX) <= opt.halfX + 1e-12);
if isempty(xg), error('psf_axis:noLines', '横向窗口太窄，P.imageX 里取不到线'); end

% 两个时间零点修正：平面的逐角度 dt、聚焦的 dt_foc（E0 标定，见 stage0_calib）
cf = fullfile(P.dataDir, 'calib.mat');
if exist(cf, 'file')
    C = load(cf);
    dt = C.dt(:);
    if isfield(C, 'dt_foc'), dtf = C.dt_foc; else, dtf = 0; end
else
    dt = zeros(P.nAngles, 1);  dtf = 0;
end
ia0 = P.zeroAngleIndex;                           % 0° 那次发射的下标（= 6）

[Th, ThR] = f2_setup(P);                          % 与正式采集相同的两个句柄

M = struct('zg', zg, 'xg', xg, 'opt', opt);
M.wPW = nan(nz,1);  M.wC11 = nan(nz,1);  M.wFoc = nan(nz,1);
M.xPW = nan(nz,1);  M.xC11 = nan(nz,1);  M.xFoc = nan(nz,1);
M.pPW = nan(nz,1);  M.pC11 = nan(nz,1);  M.pFoc = nan(nz,1);

for iz = 1:nz
    z0 = zg(iz);
    % 该深度的局部成像窗口：横向 ±5 mm、轴向 ±3 mm，与 metrics_resolution 的
    % 搜索窗口同尺寸——窗口边缘就是它能找到交点的边界。
    zs = P.imageZ(abs(P.imageZ - z0) <= opt.halfZ + 1e-12);
    pos = [0 0 z0];                               % 只有这一个散射子
    amp = opt.amp;

    % ---- PW：一次 0° 平面波发射 ---------------------------------------
    if doPW
        [rf, t0] = acquire_pw(Th, ThR, P, P.angles(ia0), pos, amp);
        A = das_pw(rf, t0, P.angles(ia0), P, xg, zs, dt(ia0));
        m = metrics_resolution(abs(hilbert(A)), xg, zs, 0, z0);
        M.wPW(iz) = m.wLat;  M.xPW(iz) = m.xpk;  M.pPW(iz) = m.peak;
    end

    % ---- C11：11 个角度发射，射频域相加后取包络 --------------------------
    % 顺序不能反：先把 11 个角度的复数射频相加，再取模，才是相干复合。
    if doC11
        A3 = zeros(numel(zs), numel(xg), P.nAngles);
        for ia = 1:P.nAngles
            [rf, t0] = acquire_pw(Th, ThR, P, P.angles(ia), pos, amp);
            A3(:,:,ia) = das_pw(rf, t0, P.angles(ia), P, xg, zs, dt(ia));
        end
        m = metrics_resolution(abs(sum(hilbert(A3), 3)), xg, zs, 0, z0);
        M.wC11(iz) = m.wLat;  M.xC11(iz) = m.xpk;  M.pC11(iz) = m.peak;
    end

    % ---- 聚焦：逐线发射，每条线只重建自己那一列 -------------------------
    if doFoc
        img = zeros(numel(zs), numel(xg));
        for k = 1:numel(xg)
            [rf, t0] = acquire_focused(Th, ThR, P, xg(k), pos, amp);
            img(:,k) = das_focused(rf, t0, xg(k), P, zs, dtf);
        end
        m = metrics_resolution(abs(hilbert(img)), xg, zs, 0, z0);
        M.wFoc(iz) = m.wLat;  M.xFoc(iz) = m.xpk;  M.pFoc(iz) = m.peak;
    end

    fprintf('    z = %4.1f mm : PW %s | C11 %s | 聚焦 %s mm\n', z0*1e3, ...
            nanstr(M.wPW(iz)), nanstr(M.wC11(iz)), nanstr(M.wFoc(iz)));
end
end

% =========================================================================
function s = nanstr(v)
%NANSTR  这一档没算（或没测到）时打印 "--" 而不是 NaN，进度表才读得下去。
if isnan(v), s = ' -- '; else, s = sprintf('%.3f', v*1e3); end
end
