function [img, w2] = das_focused(rf, t0, xf, P, z, dt)
%DAS_FOCUSED  单条聚焦 B 模式线的延时叠加波束形成。
%
%   [IMG, W2] = DAS_FOCUSED(RF, T0, XF, P, Z, DT)
%
%   输入
%     rf  ：[Nt x Nch] 一次聚焦发射的通道数据（acquire_focused 的输出）
%     t0  ：第 1 个采样点的绝对时刻 [s]
%     xf  ：该条线的横向位置 [m]
%     P   ：参数结构体
%     z   ：轴向像素坐标 [m]
%     dt  ：全局时间零点校正（E0 得到的 dt_foc），默认 0
%
%   输出
%     img ：[numel(z) x 1] 该条线的波束形成结果
%     w2  ：[numel(z) x 1] 插值权重平方和（噪声增益）
%
%   只重建落在该线上的像素（标准的逐线扫描），所以输出是列向量；这一点与
%   das_pw 不同——后者对一次平面波发射直接输出整幅二维图像。
%
%   延时模型（设计文档式 7）
%   ------------------------
%     tau = (z + sqrt(z^2 + (xf-x1)^2))/c     发射（按直线路径近似）+ 接收
%         + (|F| - |F-x1|)/c                  Field II 时间轴修正（报告 6.1）
%         + dt                                全局时间零点（E0）
%
%   其中 (|F| - |F-x1|)/c 正是 acquire_focused 施加的发射延时，原因与 das_pw
%   相同：Field II 会把它加到每个通道的时间轴上。发射部分的 z/c 是聚焦线成像
%   常用的近似（脉冲同时到达焦点，之后按声速继续传播），在焦区附近足够准确，
%   而这个参照方法本来就是在焦区做对比的。

if nargin < 6 || isempty(dt), dt = 0; end

z = z(:);                       % 强制深度为列向量
[Nt, Nch] = size(rf);
Nz = numel(z);
dxe = abs(P.elementX - xf);           % 1 x N，阵元到该线的距离
hw  = z / (2*P.fNumber);           % Nz x 1，动态接收半孔径
m   = dxe <= hw;                % Nz x N 接收孔径掩膜
off = (0:Nch-1) * Nt;           % 各通道线性索引偏移

tau = (z + sqrt(z.^2 + dxe.^2))/P.soundSpeed + dt ...
      + (sqrt(xf^2 + P.txFocusDepth^2) - sqrt((xf - P.elementX).^2 + P.txFocusDepth^2))/P.soundSpeed;
pos = (tau - t0) * P.sampleRate + 1;    % 小数采样下标
i0  = floor(pos);
fr  = pos - i0;
m   = m & (i0 >= 1) & (i0 < Nt);

i0 = min(max(i0, 1), Nt-1);     % 钳位，保证取值合法
v0 = rf(i0 + off);
v1 = rf(i0 + 1 + off);
w0 = (1 - fr) .* m;             % 线性插值权重（乘孔径掩膜）
w1 = fr .* m;

img = sum(w0.*v0 + w1.*v1, 2);  % 沿接收孔径相干求和
w2  = sum(w0.^2 + w1.^2, 2);    % 噪声功率增益
end
