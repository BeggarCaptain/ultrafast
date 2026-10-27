function P = params(root)
%PARAMS  全部仿真参数集中在此（唯一真源）。
%
%   P = PARAMS(ROOT) 返回本项目所有脚本/函数共用的参数结构体。ROOT 是工程根
%   目录（包含 config/、src/、run/、data/、figures/ 的那一层）；省略时用当前
%   工作目录。
%
%   命名约定
%   --------
%   字段名一律写成能自解释的英文，不用 N、c、fs、L 这类单字母缩写，读代码时不
%   必回头查表。三条习惯：
%     * n 开头 = 个数（nElements、nImageLines、nScatterers）；
%     * X / Z 结尾 = 坐标、位置或深度向量，单位米（elementX、imageZ、cystCenterX）；
%     * 单位写在同一行注释里，不靠字段名暗示。
%   公式符号与字段的对应：c -> soundSpeed，fs -> sampleRate，f0 -> centerFreq，
%   lambda -> wavelength，N -> nElements，F# -> fNumber，alpha -> angles。
%
%   为什么所有参数放在一个文件里？
%   ------------------------------
%   采集、波束形成、声场、指标、绘图五个阶段必须对“网格、角度、体模、脉冲”
%   有完全一致的理解。参数集中在这里，改一个数（例如接收 F 数）会同时作用到
%   所有阶段，报告里也只需要引用一张参数表。
%
%   单位约定：除字段名特别说明外（mm、dB 等），一律使用国际单位制
%   （米、秒、赫兹、弧度）；只有在画图时才换算成毫米。
%
%   参考文献：G. Montaldo, M. Tanter, J. Bercoff, N. Benech, M. Fink,
%   "Coherent Plane-Wave Compounding for Very High Frame Rate Ultrasonography
%   and Transient Elastography", IEEE TUFFC 56(3):489-506, 2009.

% -------------------------------------------------------------------------
% 工程目录结构。所有路径都由 ROOT 推导，不写死绝对路径，整个文件夹可以整体
% 拷贝到另一台机器上直接运行。
% -------------------------------------------------------------------------
if nargin < 1 || isempty(root)
    root = pwd;
end

P.rootDir   = root;
P.dataDir   = fullfile(root, 'data');
P.figureDir = fullfile(root, 'figures');
P.docsDir   = fullfile(root, 'docs');

% ---------------------------------------------------------------- 声学参数
% soundSpeed ：介质声速（论文取值，也是 Field II 默认值）
% sampleRate ：仿真采样率。4.5 MHz 下 30 MHz 约合每周期 6.7 个采样点，线性插值
%              误差可以压到 -20 dB 以下。
% centerFreq ：探头标称中心频率（论文），只用于计算波长 wavelength。
% bandwidthGoal：两程脉冲的“半功率（-3 dB）”相对带宽设计目标。论文没有给出带宽
%              指标，这是本实验自定的脉冲设计目标；-6 dB 口径的实测值只报告、
%              不设目标。
P.soundSpeed    = 1540;       % m/s
P.sampleRate    = 30e6;       % Hz，仿真采样率
P.centerFreq    = 4.5e6;      % Hz，标称中心频率
P.bandwidthGoal = 0.60;       % 半功率（-3 dB）相对带宽，设计目标（非文献值）

% ------------------------------------------------------------------- 阵列
% 128 阵元线阵（论文实际是 124 个有效阵元，差异见报告 3.4 节）。阵元宽度、
% 高度、kerf 论文没有给，按 0.33 mm 间距探头的典型值取。
%
% 仰角方向不做聚焦：本仿真里高度方向只有一个子阵元（no_sub_y = 1），无法表达
% 仰角延时剖面，xdc_linear_array 的焦点参数因此传 [0 0 0]。实测传 30 mm 与传 0
% 的 calc_scat_multi 输出逐位相同（差值为 0），所以配置里不留这个参数。
P.nElements = 128;            % 阵元数
P.elemPitch = 0.33e-3;        % m，阵元间距
P.elemWidth = 0.30e-3;        % m，阵元宽度（kerf = 间距 - 宽度 = 0.03 mm）
P.elemHeight = 5.0e-3;        % m，阵元高度（仰角方向孔径）
P.elemKerf  = 0.03e-3;        % m，阵元间隙
P.fNumber   = 1.75;           % 接收 F 数（论文），不做变迹
P.elementX  = ((0:P.nElements-1) - (P.nElements-1)/2) * P.elemPitch;   % 阵元中心坐标 [m]

% ------------------------------------------------------------ 发射脉冲
% 发射脉冲是高斯包络正弦：载频 4.50 MHz（= 论文标称频率）、sigma = 0.32 个载波
% 周期、窗长 5 个周期。波形由 src/pulse_waveform.m 统一给出，f2_setup 与图 F2
% 都调用它，避免“图里的脉冲”和“仿真用的脉冲”各写一遍而走样。
% 分数带宽只由无量纲量 sigma*fc 决定（半功率 BW = 0.1874/(sigma*fc)），0.32 是
% 60 % 半功率目标的反解值 0.312；与载频的绝对值无关。
% 两程回波在频域上是单程谱的平方再乘孔径的空间响应，中心因此下移到约 4.12 MHz
% （低 8 %）、半功率带宽约 62 %。这个下移来自“裸阵元没有电-机带通整形”，按偏差
% 记录在报告 3.1 节，不再用抬高载频的办法把它补回来。
% 论文没有给出带宽指标，60 % 是本实验按半功率口径设定的设计目标。
P.pulseSigmaCycles = 0.32;    % 高斯 sigma，单位是 1/载频（0.32 个载波周期）
P.pulseCarrierFreq = 4.50e6;  % Hz，发射脉冲载频（= 论文标称中心频率）
% 注意：脉冲数组的长度不在这里配置。它是派生量——峰位固定在第 3 个载波周期，
% 右侧留 6 sigma 尾巴，长度 = 3 + ceil(6*pulseSigmaCycles) 个周期，由 f2_setup
% 自动算出。早先版本曾把它做成独立参数 P.pulseHalfLen = 5：它与 sigma 会脱钩
% （改了 sigma 忘了改它就会静默截断脉冲），而且不可辨识（5/6/8 结果相同），
% 已删除。

% ---------------------------------------------------------- 成像网格
% 复合/显示图像是 256 线、横向 ±21.8 mm、深度 5–60 mm 的规则网格。横向取
% lambda/2，是 4.5 MHz 探头常用的奈奎斯特采样。
P.wavelength  = P.soundSpeed / P.centerFreq;   % m，标称频率下的波长
P.arrayLength = P.nElements * P.elemPitch;     % m，阵列长度（42.24 mm）
P.nImageLines = 256;                           % 图像线数（复合/显示用）
P.pixelPitchX = P.wavelength / 2;              % m，横向像素（0.171 mm）
P.depthMin    = 5e-3;                          % m，最浅成像深度
P.depthMax    = 60e-3;                         % m，最深成像深度
P.pixelPitchZ = 0.05e-3;                       % m，轴向像素（0.05 mm）
P.imageX      = ((0:P.nImageLines-1) - (P.nImageLines-1)/2) * P.pixelPitchX;  % m，各线横向坐标
P.imageZ      = (P.depthMin:P.pixelPitchZ:P.depthMax)';                     % m，深度向量

% ------------------------------------------------------- 常规聚焦 B 模式
% 每条线发射一次，焦点 30 mm，发射孔径 52 阵元（17.14 mm，发射 F 数 1.75）。
% 点目标体模用 256 线，使其 PSF 与复合图像采样一样细。
P.txFocusDepth = 30e-3;                        % m，发射焦点深度
P.nTxLines     = 128;                          % 散斑/囊肿体模的线数
P.txLinePitch  = 0.33e-3;                      % m，线间距（一个阵元一条线）
P.txLineX      = ((0:P.nTxLines-1) - (P.nTxLines-1)/2) * P.txLinePitch;      % 各线坐标
P.txApodWidth  = P.txFocusDepth / P.fNumber;   % m，发射孔径宽度（17.1 mm）
P.nPsfLines    = 256;                          % 点目标体模的线数
P.psfLineX     = ((0:P.nPsfLines-1) - (P.nPsfLines-1)/2) * P.pixelPitchX;    % 各线坐标

% --------------------------------------------------------- 平面波角度
% 实验设定：11 个平面波，范围 -10° ~ +10°，等间隔 2°。
% 采集、复合、声场图与所有显示只用这一套角度，不再保留其他档位。
P.angleStepDeg   = 2;                               % 角度间隔 [deg]
P.angleIndex     = -5:5;                            % 角度序号，0 = 垂直入射
P.angles         = P.angleIndex * P.angleStepDeg * pi/180;   % rad，-10° … +10°
P.nAngles        = numel(P.angles);                 % 11
P.zeroAngleIndex = find(P.angleIndex == 0);         % 0° 那次发射的下标（= 6）
P.compoundCounts = [1 11];                          % 参与比较的方法：PW 与 C11
P.adoptedAngleCount = 11;                           % 主配置（唯一复合设置）；无代码读取
P.nOptimalRef    = P.arrayLength / (P.wavelength * P.fNumber);  % ≈70.5，仅作参考

% ---------------------------------------------------------------- 体模
% 全部体模都在 y = 0 平面（二维仿真）。P.randomSeed 固定随机数，保证每次运行、
% 每种方法看到的散射子完全相同。
P.randomSeed     = 2009;             % 随机种子（论文发表年份）
P.nScatterers    = 20000;            % P1 均匀散射子数量
P.phantomExtentX = [-22 22] * 1e-3;  % m，横向范围（略宽于显示范围）
P.phantomExtentZ = [5 60] * 1e-3;    % m，轴向范围
P.cystRadius     = 3e-3;             % m，囊肿半径（直径 6 mm）
P.cystCenterX    = [-12 0 12] * 1e-3;   % m，囊肿中心横坐标
P.cystCenterZ    = [20 35 50] * 1e-3;   % m，囊肿中心深度（3 x 3 = 9 个囊肿）
P.pointTargetX   = [0 0 0 0 0 0 -12 12] * 1e-3;      % m，点目标横坐标
P.pointTargetZ   = [15 25 30 35 45 55 30 30] * 1e-3; % m，点目标深度
P.pointTargetAmp = 30;               % 点目标幅度（够强，无噪声时 PSF 也清晰）

% ------------------------------------------------------------------- 噪声
% 电子噪声加在“通道数据”上，每个阵元、每次发射都独立。sigma_e 由“单平面波在
% 30 mm 处 20 dB SNR”标定，之后所有方法共用同一个值。
P.pwSnrTargetDB = 20;         % dB，单平面波在 z = 30 mm 的目标 SNR
P.nMonteCarlo   = 20;         % 解析噪声公式的蒙特卡洛重复次数

% --------------------------------------------------------------- 显示
P.displayRangeDB = 60;        % dB，所有 B 模式图像的显示动态范围

% ------------------------------------------------------------ 声场图
% calc_hp 用的 1 mm 网格（41 x 56 个点）。
P.fieldMapX = (-20:1:20) * 1e-3;   % m
P.fieldMapZ = (5:1:60) * 1e-3;     % m

% --------------------------------------------------------------- E7（可选）
% 可选的 GPU 实验计时协议（本机无 CUDA 设备，未执行）。预热若干次后丢弃，
% 正式重复若干次，取中位数。
P.gpuBench.warmup = 3;             % 丢弃前 N 次
P.gpuBench.repeat = 10;            % 正式重复次数
P.gpuBench.counts = [1 11];        % E7 只比较单平面波与 11 角度复合

% ------------------------------------------------------- 外部工具箱路径
% Field II 是第三方 MATLAB 工具箱；这是整个配置里唯一与具体机器相关的项。
P.fieldIIPath = 'C:\Users\HP\Desktop\~\Field_II_ver_3_30_windows';
end
