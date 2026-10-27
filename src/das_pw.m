function [img, w2] = das_pw(rf, t0, alpha, P, x, z, dt)
%DAS_PW  单次平面波发射的延时叠加（DAS）波束形成。
%
%   [IMG, W2] = DAS_PW(RF, T0, ALPHA, P, X, Z, DT)
%
%   输入
%     rf    ：[Nt x Nch] 一次发射的通道数据（由 acquire_pw 得到）
%     t0    ：rf 第 1 个采样点的绝对时刻 [s]
%     alpha ：该次发射的偏转角 [rad]
%     P     ：参数结构体（提供 soundSpeed、sampleRate、elementX、fNumber）
%     x     ：横向像素坐标 [m]（行或列都可以）
%     z     ：轴向像素坐标 [m]
%     dt    ：该角度的时间零点校正 [s]（E0 标定），默认 0
%
%   输出
%     img   ：[numel(z) x numel(x)] 动态接收聚焦后的实数射频。包络与对数压缩
%             由调用者稍后完成，这样多个角度可以在射频域相加（相干复合）。
%     w2    ：[numel(z) x numel(x)] 插值权重的平方和。线性波束形成是独立噪声
%             样本的加权求和，因此输出噪声方差 = sigma_e^2 * w2。SNR 指标就是
%             用这个关系解析传播噪声的（并在 stage_sigma 里用蒙特卡洛验证过）。
%
%   延时模型（设计文档式 1–4）
%   --------------------------
%     tau = (z*cos(alpha) + x*sin(alpha))/c   发射：平面波到达 (x,z)
%         + sqrt(z^2 + (x-x1)^2)/c            接收：回波从 (x,z) 到阵元 x1
%         + x1*sin(alpha)/c                   Field II 时间轴修正（报告 6.1）
%         + dt                                全局时间零点（E0）
%
%   接收孔径：矩形窗，|x1 - x| <= z/(2F)，不变迹（论文也不用变迹）。窗口外的
%   阵元权重为 0，这也是 w2 随深度变化的原因。

if nargin < 7 || isempty(dt), dt = 0; end

z = z(:);  x = x(:).';                  % 强制：深度为列，横向为行
[Nt, Nch] = size(rf);
Nz = numel(z);  Nx = numel(x);
img = zeros(Nz, Nx);
w2  = zeros(Nz, Nx);

ca = cos(alpha);  sa = sin(alpha);
off = (0:Nch-1) * Nt;                   % 各通道在线性索引中的偏移
dxe = P.elementX;                       % 1 x N 阵元坐标

for k = 1:Nx                            % 一次算一条线（内存受限）
    xc = x(k);
    axc = abs(dxe - xc);                % 阵元到像素的横向距离，1 x N
    hw  = z / (2*P.fNumber);            % 各深度的半孔径宽度，Nz x 1
    m   = axc <= hw;                    % Nz x N 孔径掩膜（隐式扩展）

    % Field II 返回的每个通道都把自己的发射延时加在了时间轴上（E0 实测：
    % 带此项残差 0.01 µs，不带则 1.2 µs）。所以延时必须包含这一项。
    tau = (z*ca + xc*sa)/P.soundSpeed + sqrt(z.^2 + (dxe - xc).^2)/P.soundSpeed ...
          + dxe*sa/P.soundSpeed + dt;
    pos = (tau - t0) * P.sampleRate + 1;   % 1 基的（小数）采样下标

    i0 = floor(pos);                    % 左侧相邻采样
    fr = pos - i0;                      % 插值权重，[0,1)
    m  = m & (i0 >= 1) & (i0 < Nt);     % 延时超出记录范围的像素置零

    i0 = min(max(i0, 1), Nt-1);         % 钳位，保证取值下标永远合法
    v0 = rf(i0 + off);                  % 左侧采样值
    v1 = rf(i0 + 1 + off);              % 右侧采样值

    w0 = (1 - fr) .* m;                 % 线性插值权重（已乘掩膜）
    w1 = fr .* m;

    img(:, k) = sum(w0.*v0 + w1.*v1, 2);% 沿阵元方向相干求和
    w2(:, k)  = sum(w0.^2 + w1.^2, 2);  % 同一求和的噪声功率增益
end
end
