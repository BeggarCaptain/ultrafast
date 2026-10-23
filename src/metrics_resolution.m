function m = metrics_resolution(env, x, z, xp, zp)
%METRICS_RESOLUTION  单个点目标的 -10 dB 横向宽度与 -6 dB 轴向宽度。
%
%   M = METRICS_RESOLUTION(ENV, X, Z, XP, ZP) 测量标称位置在 (XP, ZP) 的目标
%   的点扩散函数。
%
%   输入
%     env ：包络图像 [numel(z) x numel(x)]
%     x   ：横向坐标向量 [m]
%     z   ：轴向坐标向量 [m]
%     xp  ：目标的标称横向位置 [m]
%     zp  ：目标的标称轴向位置 [m]
%
%   输出结构体
%     w10  ：横向 -10 dB 全宽 [m]
%     w6   ：轴向 -6 dB 全宽 [m]
%     xpk  ：实测峰值横向位置 [m]
%     zpk  ：实测峰值轴向位置 [m]
%     peak ：峰值处的包络值（用于合理性检查）
%
%   步骤
%   ----
%   1. 在标称位置附近取 ±5 mm / ±3 mm 的窗口（足够容纳整个主瓣，又不至于把
%      相邻目标框进来）；
%   2. 求窗口内的全局最大值；
%   3. 取峰值所在行的横向曲线、所在列的轴向曲线；
%   4. 两条曲线都用 pchip 上采样 8 倍，再用线性插值定位 -10 dB / -6 dB 交叉点。
%
%   第 4 步很关键：轴向 PSF 只有约 0.18 mm 宽，而图像网格是 0.05 mm，-6 dB
%   以上只有三个采样点左右。不上采样的话宽度会被量化成 0.05 mm 的整数倍，
%   在不同深度、不同方法之间来回跳。用 pchip 而不是 spline，是因为 pchip 不会
%   在峰值处过冲（过冲会把测得的宽度压小）。

ix = find(abs(x - xp) <= 5e-3);         % 目标附近的横向窗口
iz = find(abs(z - zp) <= 3e-3);         % 目标附近的轴向窗口
E  = double(env(iz, ix));
[~, im] = max(E(:));                    % 窗口内的峰值
[r, c] = ind2sub(size(E), im);

% ---- 过峰值、上采样后的两条曲线 ----------------------------------------
xu = linspace(x(ix(1)), x(ix(end)), 8*numel(ix));
zu = linspace(z(iz(1)), z(iz(end)), 8*numel(iz));
lat = interp1(x(ix).', E(r, :).', xu, 'pchip');   % 横向曲线（第 r 行）
axp = interp1(z(iz),    E(:, c),    zu, 'pchip'); % 轴向曲线（第 c 列）

% ---- 峰值位置与宽度 -----------------------------------------------------
[~, iL] = max(lat);  m.xpk = xu(iL);
[~, iA] = max(axp);  m.zpk = zu(iA);
m.w10 = crossing_width(xu, lat, 10^(-10/20));     % 横向 -10 dB 宽度
m.w6  = crossing_width(zu, axp, 0.5);             % 轴向 -6 dB 宽度
m.peak = max(lat);
end
