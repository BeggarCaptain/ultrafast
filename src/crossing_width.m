function w = crossing_width(v, p, lvl)
%CROSSING_WIDTH  一维曲线在指定电平处的宽度（交叉点线性插值）。
%
%   W = CROSSING_WIDTH(V, P, LVL) 返回曲线 P 在峰值两侧穿过 LVL*max(P) 的两点
%   之间的距离。交叉点先在采样网格上定位，再在相邻两个采样点之间做线性插值
%   细化，因此结果不会被像素尺寸量化（这一点很关键：轴向 PSF 只有约 0.19 mm
%   宽，而网格是 0.05 mm）。
%
%   输入
%     v   ：坐标向量（mm 或 m，返回值单位与之一致）
%     p   ：曲线（包络或强度），长度与 V 相同
%     lvl ：相对峰值的线性阈值，例如
%             10^(-10/20) = 0.316  求 -10 dB 宽度
%             0.5                  求 -6 dB 宽度
%
%   输出
%     w   ：该电平处的全宽；若某一侧始终没降到阈值以下（例如平顶曲线），返回 NaN。

v = v(:);
p = p(:);

[~, ip] = max(p);                       % 峰值位置
thr = lvl * p(ip);                      % 换算成绝对阈值
n = numel(p);

% ---- 从峰值向左走，直到曲线降到阈值以下 --------------------------------
iL = ip;
while iL > 1 && p(iL) > thr
    iL = iL - 1;
end
if iL == ip                            % 峰值本身已在阈值以下
    w = NaN; return;
end
vL = interp1([p(iL); p(iL+1)], [v(iL); v(iL+1)], thr, 'linear', 'extrap');

% ---- 从峰值向右走，同理 ------------------------------------------------
iR = ip;
while iR < n && p(iR) > thr
    iR = iR + 1;
end
if iR == ip                            % 右半边没有穿过阈值
    w = NaN; return;
end
vR = interp1([p(iR-1); p(iR)], [v(iR-1); v(iR)], thr, 'linear', 'extrap');

w = vR - vL;                           % 该电平处的全宽
end
