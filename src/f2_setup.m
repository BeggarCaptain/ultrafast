function Th = f2_setup(P)
%F2_SETUP  初始化 Field II 并建立 128 阵元线阵。
%
%   Th = F2_SETUP(P) 返回换能器句柄，激励与冲激响应都已装载，可直接调用
%   calc_scat_multi / calc_hp。
%
%   本函数可以在同一个 MATLAB 会话里反复调用：每次都会重新初始化 Field II
%   并从头建立阵列，因此采集、噪声、声场各阶段都从确定的状态开始。
%
%   重要：这里刻意从不调用 field_end。在本机这套 Field II 版本上，
%   “field_end 之后再 field_init”会把 MATLAB 进程直接搞崩（曾经杀死过交互会话）；
%   单独调用 field_init(-1) 已经足以（重新）初始化，所以只用它。
%
%   脉冲设计
%   --------
%   波束形成看到的“两程”脉冲，近似等于换能器冲激响应被作用两次，所以单程高斯
%   必须比目标结果更短、载频更高。取 sigma = 0.32/f0、载频 4.90 MHz 时，两程
%   脉冲实测为 4.45 MHz / 62.8 % 带宽，与论文一致（E0 会验证）。

% ---- 确保 Field II 工具箱在搜索路径上 -----------------------------------
if exist('field_init', 'file') ~= 2
    addpath(P.fieldIIPath);
end

% ---- 库级设置 -----------------------------------------------------------
field_init(-1);             % -1 表示关闭 Field II 的图形输出
set_field('c', P.c);        % 所有延时计算使用的声速
set_field('fs', P.fs);      % 返回信号的采样率
set_field('show_times', 0); % 不打印各步骤耗时

% ---- 线阵几何 -----------------------------------------------------------
% 参数依次为：阵元数、宽度、高度、kerf、x 方向子阵元数、y 方向子阵元数、焦点。
% 子阵元数取 1，表示每个物理阵元只用一个数学子阵元：对 0.3 mm 阵元、从 5 mm
% 开始成像的场景，这是速度与精度之间较好的折中。
%
% 第七个参数是固定的仰角焦点，这里传 [0 0 0]（不聚焦）。原因：y 方向只有 1 个
% 子阵元，高度方向无法产生延时剖面；实测把该值改成 30 mm 后，三个不同位置
% （轴上与离轴）的 calc_scat_multi 输出与原配置逐位相同，最大差值为 0。既然
% 它从不起作用，配置里就不保留这个参数。
Th = xdc_linear_array(P.N, P.width, P.height, P.kerf, 1, 1, [0 0 0]);

% ---- 脉冲：高斯包络正弦，存为冲激响应 -----------------------------------
fc  = P.pulseCarrier;                              % 单程脉冲载频 [Hz]
t   = (0:round(P.pulseHalfLen/fc*P.fs))' / P.fs;   % 脉冲时间轴 [s]
sig = P.pulseSigma / fc;                           % 高斯 sigma [s]
imp = exp(-((t - 3/fc).^2) / (2*sig^2)) .* sin(2*pi*fc*(t - 3/fc));

xdc_impulse(Th, imp');      % Field II 要求行向量，这里转置
xdc_excitation(Th, 1);      % 单位冲激激励：单程脉冲就等于 imp

% ---- 全孔径、不变迹（与论文一致），焦点归零 -----------------------------
xdc_apodization(Th, 1, ones(1, P.N));   % 发射：所有阵元参与
xdc_apodization(Th, 2, ones(1, P.N));   % 接收：硬件不做加窗
xdc_center_focus(Th, [0 0 0]);          % 中性初始焦点
end
