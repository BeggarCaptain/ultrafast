function [rf, t0] = acquire_pw(Th, ThR, P, alpha, pos, amp)
%ACQUIRE_PW  一次平面波发射，返回 128 通道原始射频（波束形成之前）。
%
%   [RF, T0] = ACQUIRE_PW(TH, THR, P, ALPHA, POS, AMP)
%
%   输入
%     Th    ：发射句柄（由 f2_setup 生成，可加窗、可聚焦）
%     ThR   ：接收句柄（同一几何，变迹恒为单位）——平面波本来就全孔径，两个
%             句柄的设置此时完全一样，分开只是为了和聚焦模式共用一条代码路径
%     P     ：参数结构体（config/params.m）
%     alpha ：偏转角 [rad]，取自论文式 (5) 的角度集合
%     pos   ：N x 3 散射子坐标 [m]（二维体模 y = 0）
%     amp   ：N x 1 散射子幅度
%
%   输出
%     rf    ：[Nt x 128] 每个阵元的接收信号（即自写 DAS 要处理的“通道数据”）
%     t0    ：第一个采样点对应的绝对时刻 [s]；第 k 个采样点对应 t0 + k/fs
%
%   为什么一次平面波发射就能覆盖整幅图
%   ----------------------------------
%   平面波没有发射焦点：一次发射同时照亮整个视场，所以一帧只需 1 次发射（聚焦
%   B 模式要 128 次）。代价是发射端不聚焦，横向分辨率完全靠接收端逐像素补偿
%   （das_pw 的延时求和），因此单平面波的分辨率与对比度都低于聚焦发射——这正是
%   论文要用多角度复合去补的那部分。
%
%   平面波是这样产生的：每个阵元按 d(x1) = x1*sin(alpha)/c 的线性延时发射，
%   相邻阵元的延时差对应波前沿 x 方向倾斜 alpha（等相位面是斜直线，不是圆弧）。
%   符号在 E0 中验证过：实测峰值时刻随 x 的斜率与 +sin(alpha)/c 相差不到 2 %。
%
%   给自写 DAS 的提醒：Field II 会把 d(x1) 加到第 x1 通道的时间轴上（实测：
%   给接收句柄单独设延时，该通道波形整体平移 d(x1)），所以波束形成的延时里
%   必须包含同一项（+x1*sin(alpha)/c），否则图像会横向剪切——见 das_pw 与
%   docs/实验报告_讨论与附录.md 的 6.1 节。

d = P.elementX * sin(alpha) / P.soundSpeed;              % 各阵元发射延时 [s]，可能为负

% 两个句柄的延时必须设成同一个 d：Field II 把它写进每条通道的时间轴，DAS 里
% 再减掉，两边不一致图像就错位。参数 0 是“从 t=0 起生效”的时间线 [s]，不是
% 发射/接收开关（见 f2_setup 的说明）。
xdc_center_focus(Th, [0 0 0]);            % 清掉上一次的焦点
xdc_focus_times(Th, 0, d);                % 发射延时
xdc_apodization(Th, 0, ones(1, P.nElements));     % 全孔径、不变迹（论文设置）

xdc_center_focus(ThR, [0 0 0]);
xdc_focus_times(ThR, 0, d);               % 与发射端相同的时间基准
xdc_apodization(ThR, 0, ones(1, P.nElements));    % 接收永远不变迹

[rf, t0] = calc_scat_multi(Th, ThR, pos, amp);   % 脉冲回波，每个阵元一列
end
