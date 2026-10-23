function [rf, t0, d] = acquire_focused(Th, P, xf, pos, amp)
%ACQUIRE_FOCUSED  一次聚焦发射，返回 128 通道原始射频。
%
%   [RF, T0, D] = ACQUIRE_FOCUSED(TH, P, XF, POS, AMP)
%
%   输入
%     Th   ：Field II 换能器句柄（由 f2_setup 生成）
%     P    ：参数结构体
%     xf   ：该条扫描线的横向位置 [m]
%     pos  ：N x 3 散射子坐标 [m]
%     amp  ：N x 1 散射子幅度
%
%   输出
%     rf   ：[Nt x 128] 通道数据
%     t0   ：第一个采样点的绝对时刻 [s]
%     d    ：本次实际使用的发射延时向量 [s]（行向量）
%
%   发射孔径
%   --------
%   以 XF 为中心、宽 P.txApodW = 17.14 mm 的矩形窗，即 52 个阵元，焦点处
%   F 数 = 30/17.14 = 1.75（论文取值）。靠近阵列两端时窗口会被物理孔径截断，
%   这就是设计里说的“边缘线自动裁剪”。
%
%   为什么这里要显式给延时
%   ----------------------
%   xdc_focus 也能算出同样的延时，但它不返回这个向量。DAS 需要一模一样的
%   向量，因为 Field II 会把发射延时加到每个通道的时间轴上（E0 实测：带这一项
%   残差 0.01 µs，不带则 1.2 µs）。在这里算出 d、再在 das_focused 里加回去，
%   发射与接收的时间基准才一致。

F = [xf 0 P.txFocus];                                  % 焦点坐标 [m]
% d(x1) = (|F| - |F - x1|)/c：离焦点越远的阵元越早发射
d = (norm(F) - sqrt((P.xe - xf).^2 + P.txFocus^2)) / P.c;

xdc_center_focus(Th, [0 0 0]);                         % 清掉上一次的焦点
xdc_focus_times(Th, 1, d);                             % 1：仅发射

apo = double(abs(P.xe - xf) <= P.txApodW/2 + 1e-12);   % 52 阵元矩形窗，边缘自动截断
xdc_apodization(Th, 1, apo);                           % 施加在发射上

[rf, t0] = calc_scat_multi(Th, Th, pos, amp);          % 脉冲回波通道数据
end
