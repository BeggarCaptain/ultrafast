function s = pulse_stats(rf, fs)
%PULSE_STATS  测量两程脉冲的中心频率与 -6 dB 带宽。
%
%   S = PULSE_STATS(RF, FS) 从“单散射子”采集的数据测量脉冲回波小波的频谱：
%   取最强通道、加窗、补零、做 FFT，然后定位 -6 dB 的两个边沿。
%
%   输入
%     rf ：[Nt x Nch] 单个点散射子的通道数据
%     fs ：采样率 [Hz]
%
%   输出结构体
%     f1, f2 ：下/上 -6 dB 频率 [Hz]
%     fc     ：中心频率，(f1+f2)/2 [Hz]
%     bw     ：相对带宽，(f2-f1)/fc
%     f, S   ：频率轴 [Hz] 与频谱 [dB]，供图 F2 使用
%     wav    ：实际参与分析的加窗波形
%     ich    ：所用的通道号
%
%   为什么取最强通道：点散射子在阵列轴线上，各通道看到的脉冲基本相同，取最强
%   的那个只是为了让测量信噪比最高。记录很短（约 300 个采样点），所以用 Tukey
%   窗抑制截断旁瓣，并用 8 倍补零把 -6 dB 交叉点插值出来。
%
%   只保留 0..fs/2 这一半：实信号的 FFT 是对称的，如果整轴搜索，会先碰到
%   fs - f 处的镜像，报出一个荒唐的带宽（这个坑真的踩过一次）。

[~, ich] = max(max(abs(rf), [], 1));        % 最强通道
w    = rf(:, ich) .* tukeywin(size(rf, 1), 0.1);   % 端点加窗，抑制截断旁瓣
nfft = 8 * numel(w);                        % 8 倍补零，得到更细的频率网格
S    = abs(fft(w, nfft));
f    = (0:nfft-1)' * fs / nfft;

half = f < fs/2;                            % 只保留物理上有意义的一半
f = f(half);  S = S(half);
S = S / max(S);                             % 归一化到 0 dB

idx = find(S >= 10^(-3/20));                % -6 dB 点（幅度比 0.708）
s.f1 = f(idx(1));                           % 下边沿
s.f2 = f(idx(end));                         % 上边沿
s.fc = (s.f1 + s.f2) / 2;                   % 中心频率
s.bw = (s.f2 - s.f1) / s.fc;                % 相对 -6 dB 带宽
s.f  = f;
s.S  = 20*log10(max(S, 1e-12));             % dB 频谱，加下限以便画图
s.wav = w;
s.ich = ich;
end
