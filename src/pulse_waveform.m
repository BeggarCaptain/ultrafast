function [w, t] = pulse_waveform(P)
%PULSE_WAVEFORM  发射脉冲（单程）的时域波形：高斯包络正弦。
%
%   [W, T] = PULSE_WAVEFORM(P) 返回列向量 W 与时间轴 T [s]，采样率为
%   P.sampleRate。W 的峰值归一化前的最大绝对值约为 1（正弦幅度为 1）。
%
%   这是全项目唯一的脉冲定义。f2_setup 把它交给 xdc_impulse，图 F2 直接画它，
%   两处共用同一份代码，改载频或改 sigma 时不会出现“图与仿真不是同一个脉冲”。
%
%   两个约定（都是工程约定，不是推导结果）
%   --------------------------------------
%   峰位固定在第 3 个载波周期：沿袭早期写法，实际只要头部余量 ≥5σ≈1.6 个周期
%   就够，取 3 只是为了和旧结果可比。右侧留 6σ 的尾巴（残余约 1e-8，5σ 也够）。
%   窗长是派生量：3 + ceil(6σ) 个载波周期，不可辨识——5/6/8 个周期给出同一条
%   脉冲。早先它曾是独立参数
%   P.pulseHalfLen，会和 sigma 脱钩（改了 sigma 忘了改它就会截断脉冲），已删除。
%
%   为什么峰值不放在第 1 个周期：头部余量只有 1/σ=3.1σ 时 t=0 处的包络为 0.76 %
%   （σ=0.32 时还能接受），但把 σ 调大到 0.5 就会变成 13.5 %，截断开始影响带宽
%   读数（实测 −6 dB 带宽偏高 6.6 个百分点）。3 个周期在 σ ≤ 0.5 时都保持 ≥6σ。

fc  = P.pulseCarrierFreq;                        % 载频 [Hz]
sig = P.pulseSigmaCycles / fc;                   % 高斯 sigma [s]
nCyc = 3 + ceil(6*P.pulseSigmaCycles);           % 窗长 [载波周期数]

t = (0:round(nCyc/fc*P.sampleRate))' / P.sampleRate;          % 时间轴 [s]
w = exp(-((t - 3/fc).^2) / (2*sig^2)) .* sin(2*pi*fc*(t - 3/fc));
end
