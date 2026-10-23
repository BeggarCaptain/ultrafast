function [rf, t0] = acquire_pw(Th, P, alpha, pos, amp)
%ACQUIRE_PW  一次平面波发射，返回 128 通道原始射频（波束形成之前）。
%
%   [RF, T0] = ACQUIRE_PW(TH, P, ALPHA, POS, AMP)
%
%   输入
%     Th    ：Field II 换能器句柄（由 f2_setup 生成）
%     P     ：参数结构体（config/params.m）
%     alpha ：偏转角 [rad]，取自论文式 (5) 的角度集合
%     pos   ：N x 3 散射子坐标 [m]（二维体模 y = 0）
%     amp   ：N x 1 散射子幅度
%
%   输出
%     rf    ：[Nt x 128] 每个阵元的接收信号（即自写 DAS 要处理的“通道数据”）
%     t0    ：第一个采样点对应的绝对时刻 [s]；第 k 个采样点对应 t0 + k/fs
%
%   平面波是这样产生的：每个阵元按 d(x1) = x1*sin(alpha)/c 的线性延时发射，
%   波前因此倾斜 alpha。符号在 E0 中验证过：实测峰值时刻随 x 的斜率与
%   +sin(alpha)/c 相差不到 2 %。
%
%   给自写 DAS 的提醒：Field II 会把 d(x1) 加到第 x1 通道的时间轴上，所以
%   波束形成的延时里必须包含同一项（+x1*sin(alpha)/c），否则图像会横向剪切
%   ——见 das_pw 与报告 6.1 节。

d = P.xe * sin(alpha) / P.c;              % 各阵元发射延时 [s]，可能为负

xdc_center_focus(Th, [0 0 0]);            % 清掉上一次的焦点
xdc_focus_times(Th, 0, d);                % 0：同时设发射（与接收）延时
xdc_apodization(Th, 1, ones(1, P.N));     % 全孔径、不变迹（论文设置）

[rf, t0] = calc_scat_multi(Th, Th, pos, amp);   % 脉冲回波，每个阵元一列
end
