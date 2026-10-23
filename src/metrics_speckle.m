function s = metrics_speckle(env, x, z, xc, zc, half)
%METRICS_SPECKLE  散斑一阶统计量与相关长度。
%
%   S = METRICS_SPECKLE(ENV, X, Z, XC, ZC, HALF)
%
%   输入
%     env  ：包络图像 [numel(z) x numel(x)]
%     x,z  ：坐标向量 [m]
%     xc,zc：分析 ROI 的中心 [m]
%     half ：方形 ROI 的半边长 [m]，默认 5 mm（即 10 x 10 mm）
%
%   输出结构体
%     ratio  ：强度的 μ/σ（理想指数散斑 -> 1.00）
%     ratioA ：包络的 μ/σ（理想 Rayleigh 散斑 -> 1.91）
%     corrX, corrZ ：强度自相关的 -6 dB 宽度 [m]
%     hist/histEdges：归一化强度直方图，供图 F13 使用
%
%   散斑统计只有在去掉“缓慢亮度趋势”之后才有意义。仿真图像没有 TGC，10 mm 的
%   ROI 内平均强度会变化好几个 dB；这个趋势会主导方差，即使散斑完全发展，测得
%   的 μ/σ 也会掉到 0.5 左右（这是实测到的现象，也是下面要除以局部均值的原因）。
%   除以约 2 mm 的盒式平均，只保留比它细得多的散斑起伏，把趋势去掉。
%
%   关于两个比值：设计文档最初把 1.91 当作“散斑 μ/σ”的目标值。1.91 是
%   RAYLEIGH 值，对应的是包络；完全发展的散斑其“强度”服从指数分布，μ/σ = 1.00。
%   这里两个都返回，报告里各对各的理想值。

if nargin < 6 || isempty(half), half = 5e-3; end

I = double(env).^2;                                 % 强度图

% ---- 去掉缓慢的深度/横向亮度趋势 ----------------------------------------
kz = max(3, 2*round(1e-3/(z(2)-z(1))) + 1);         % 约 2 mm 的轴向窗口
kx = max(3, 2*round(1e-3/(x(2)-x(1))) + 1);         % 约 2 mm 的横向窗口
B  = imfilter(I, ones(kz,kx)/(kz*kx), 'replicate'); % 局部均值（replicate 防止边缘偏暗）
In = I ./ max(B, realmin);                          % 去趋势后的强度，均值约 1

[XX, ZZ] = meshgrid(x, z);
sel = abs(XX - xc) <= half & abs(ZZ - zc) <= half;  % 方形 ROI
v = In(sel);

s.mu     = mean(v);
s.sigma  = std(v);
s.ratio  = s.mu / s.sigma;                          % 强度：指数分布 -> 1.00
s.ratioA = mean(sqrt(v)) / std(sqrt(v));            % 包络：Rayleigh -> 1.91
s.n      = numel(v);

% ---- 相关长度：沿一个方向平均得到低方差的一维曲线，再求自相关宽度 ------
% 沿一个轴做平均后，曲线的自相关宽度就是另一个方向上的散斑尺度。
rows = any(sel, 2);  cols = any(sel, 1);
profZ = mean(In(rows, cols), 2);
profX = mean(In(rows, cols), 1);

s.corrZ = corr_len(z(rows), profZ);                 % 轴向散斑长度
s.corrX = corr_len(x(cols), profX);                 % 横向散斑长度

% ---- 归一化的一维强度直方图 --------------------------------------------
s.histEdges = linspace(0, 4, 41);                   % I/<I> 从 0 到 4（40 个区间）
s.hist = histcounts(v / s.mu, s.histEdges);
s.hist = s.hist / sum(s.hist);                      % 每个区间的概率
end

function w = corr_len(v, p)
%CORR_LEN  一维曲线自相关的 -6 dB 宽度。
%   V 是坐标向量，P 是曲线。若曲线噪声过大、自相关始终降不到 0.5 以下
%   （也就是采样网格分辨不出散斑），返回 NaN。

v = v(:);  p = p(:) - mean(p);          % 去掉直流分量
ac = xcorr(p, 'coeff');                 % 归一化自相关，共 2N-1 个延迟
ac = ac(numel(p):end);                  % 只保留延迟 0..N-1
vv = (0:numel(ac)-1)' * (v(2)-v(1));    % 延迟换算成距离

if ac(1) < 0.5
    w = NaN; return;                    % 连零延迟都低于 0.5
end
i = find(ac < 0.5, 1);                  % 第一个低于 -6 dB 的延迟
if isempty(i) || i == 1
    w = NaN; return;                    % 始终不穿越，或立刻穿越
end
w = interp1([ac(i-1); ac(i)], [vv(i-1); vv(i)], 0.5, 'linear', 'extrap');
end
