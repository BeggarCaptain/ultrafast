function make_figures()
%MAKE_FIGURES  由已保存的标定、声场与分析结果生成图 F1–F13。
%   在 run_tail / run_refresh 结尾调用。
%
%   输入（都在磁盘上，这里不重算任何仿真）
%     data/calib.mat          ：脉冲、dt、E0 指标              （stage0_calib）
%     data/field/maps.mat     ：发射声场 / 合成声场            （stage_fields）
%     data/res/results.mat    ：图像与指标                     （stage_analysis）
%
%   输出：figures/F1...F13.png，200 dpi。每张图保存后立即关闭，脚本可以无人值守
%   地跑完。
%
%   约定：所有坐标轴单位都是 mm；所有灰度 B 模式共用同一个 60 dB 色标
%   （P.drDB），这样各子图可以直接用眼睛比较。
%
%   图片清单
%   --------
%   F1  几何示意（阵列 + 三个波前）              —— 纯绘制
%   F2  两程脉冲波形与频谱                       —— calib
%   F3  单个 0° 平面波的发射声场                 —— maps
%   F4  九个偏转角的声场，显示波前倾斜           —— maps
%   F5  合成发射声场（11 个角度叠加）            —— maps
%   F6  点目标的 128 通道原始射频                —— P3_a006（0°）
%   F7  取包络前的射频：PW / C11                  —— results
%   F8  囊肿体模 B 模式，各档位 + 聚焦           —— results
%   F9  散斑体模 B 模式，各档位 + 聚焦           —— results
%   F10 分辨率随深度（横向与轴向）               —— results
%   F11 对比度与 CNR 随角度数                    —— results
%   F12 相对 SNR 随深度                          —— results
%   F13 散斑直方图、相关长度、μ/σ                —— results

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);
C = load(fullfile(P.dataDir, 'calib.mat'));
M = load(fullfile(P.dataDir, 'field', 'maps.mat'));
R = load(fullfile(P.dataDir, 'res', 'results.mat'));

F = P.figDir;
if ~exist(F, 'dir'), mkdir(F); end

% --------------------------------------------------------------------- F1
% 示意图，不是仿真结果：说明“偏转角 alpha 的平面波”在几何上是什么样子，
% 这样 F3/F4 的声场图才看得懂。
fig = figure('Position', [100 100 900 480]); hold on;

% 阵列：每四个阵元画一个小矩形（128 个全画就是一条黑杠）
for k = 1:4:P.N
    rectangle('Position', [(P.xe(k)-P.width/2)*1e3, -0.35, P.width*1e3, 0.35], ...
        'FaceColor', [0.35 0.35 0.35], 'EdgeColor', 'none');
end

% 三个波前：本实验的最大角 ±10° 与 0°
for a = [-10 0 10]*pi/180
    x0 = -25*sin(a);
    plot([x0 x0+50*sin(a)], [0 60*cos(a)], 'LineWidth', 1.2);
end
plot([0 0], [0 60], 'k--');                       % 垂直入射的参考线

text(-24, 4, 'array, 128 el., 0.33 mm pitch', 'FontSize', 9);
text(6, 55, '\alpha = +10\circ', 'FontSize', 9);
text(6, 57, '\alpha = 0\circ', 'FontSize', 9);
text(-24, 55, '\alpha = -10\circ', 'FontSize', 9);
xlabel('x (mm)'); ylabel('z (mm)'); axis([-25 25 -2 62]); grid on;
title('Plane-wave transmit geometry');
save_fig(fig, fullfile(F, 'F1_geometry.png')); close(fig);

% --------------------------------------------------------------------- F2
% 脉冲的时域与频域。-6 dB 带宽与中心频率就是 E0 的验收依据：
% 4.5 ± 0.1 MHz、60 ± 5 %。
fig = figure('Position', [100 100 900 340]);

% (a) 波形：这里的包络宽度决定了轴向分辨率（0.18 mm）
subplot(1,2,1);
tw = (0:numel(C.pulseWav)-1)'/P.fs*1e6;           % µs
plot(tw, C.pulseWav / max(abs(C.pulseWav)), 'LineWidth', 1.2);
xlabel('t (\mus)'); ylabel('normalised'); grid on;
title(sprintf('Two-way pulse, f_c = %.2f MHz', C.pulse.fc/1e6));

% (b) 频谱，并标出 -6 dB 的两个边沿
subplot(1,2,2);
plot(C.pulseF/1e6, C.pulseS, 'LineWidth', 1.2); hold on;
plot([C.pulse.f1 C.pulse.f2]/1e6, [-6 -6], 'r-', 'LineWidth', 1.5);   % -6 dB 带
plot([C.pulse.fc C.pulse.fc]/1e6, [-30 0], 'k:');                     % 中心频率
xlim([2 8]); ylim([-45 2]); grid on;
xlabel('f (MHz)'); ylabel('dB');
title(sprintf('-6 dB BW = %.1f %%', 100*C.pulse.bw));
save_fig(fig, fullfile(F, 'F2_pulse.png')); close(fig);

% --------------------------------------------------------------------- F3
% 0° 单平面波：没有焦点，近场最亮、随深度衰减（约 1/sqrt(z)）。这是与 F5 的
% 合成声场作对比的“之前”图。
i0 = find(P.iAng == 0);
mp = double(M.peakAll(:,:,i0)); mp = mp / max(mp(:));
fig = figure('Position', [100 100 950 620]);

% (a) 二维声场，-30 dB 色标（实际动态范围远大于 30 dB）
subplot(2,2,[1 3]);
imagesc(M.x*1e3, M.z*1e3, 20*log10(max(mp, 1e-4))); axis image; colormap(gca, hot);
caxis([-30 0]); colorbar; xlabel('x (mm)'); ylabel('z (mm)');
title('Transmit pressure, \alpha = 0\circ (dB)');

% (b) 30 mm 处横向切片：基本是平的，即没有横向约束
subplot(2,2,2);
[~, iz] = min(abs(M.z - 30e-3));
plot(M.x*1e3, 20*log10(max(mp(iz,:), 1e-4)), 'LineWidth', 1.2);
xlabel('x (mm)'); ylabel('dB'); grid on; title('lateral at z = 30 mm');

% (c) 轴线上的纵向切片：单调衰减，没有焦峰
subplot(2,2,4);
[~, ix] = min(abs(M.x));
plot(M.z*1e3, 20*log10(max(mp(:,ix), 1e-4)), 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('dB'); grid on; title('axial at x = 0');
save_fig(fig, fullfile(F, 'F3_field_pw0.png')); close(fig);

% --------------------------------------------------------------------- F4
% 六个偏转角（-10°/-6°/-2°/+2°/+6°/+10°）。每幅按自身最大值归一化，这样
% 波前倾斜（方向随 alpha 变号）才看得出来；这是角度实现正确的几何证据。
isel = [1 3 5 7 9 11];
fig = figure('Position', [100 100 900 500]);
for q = 1:6
    subplot(2,3,q);
    mp = double(M.peakAll(:,:,isel(q))); mp = mp / max(mp(:));
    imagesc(M.x*1e3, M.z*1e3, 20*log10(max(mp, 1e-4))); axis image;
    colormap(gca, hot); caxis([-30 0]); set(gca, 'XTick', -20:20:20);
    title(sprintf('\\alpha = %.0f\\circ', P.angles(isel(q))*180/pi), 'FontSize', 8);
    if mod(q-1,3) == 0, ylabel('z (mm)'); end
    if q >= 4, xlabel('x (mm)'); end
end
save_fig(fig, fullfile(F, 'F4_field_angles.png')); close(fig);

% --------------------------------------------------------------------- F5
% 声场部分的核心图：把 11 个角度的声场按“到焦点的传播时间”平移后相加。
% 30 mm 处出现焦斑，横向剖面收敛到物理聚焦发射的那一条——这正是“复合在
% 合成一个发射焦点”的含义。
syn = double(M.synPeak); syn = syn / max(syn(:));
foc = double(M.focPeak); foc = foc / max(foc(:));
pw1 = double(M.peakAll(:,:,i0)); pw1 = pw1 / max(pw1(:));
fig = figure('Position', [100 100 950 620]);

% (a) 二维合成声场
subplot(2,2,[1 3]);
imagesc(M.x*1e3, M.z*1e3, 20*log10(max(syn, 1e-4))); axis image; colormap(gca, hot);
caxis([-30 0]); colorbar; xlabel('x (mm)'); ylabel('z (mm)');
title(sprintf('Synthetic transmit field, %d angles (focus gain %.0f)', ...
    P.nAngles, M.gain));

% (b) 30 mm 处横向切片：合成场 vs 聚焦 vs 单平面波
subplot(2,2,2);
plot(M.x*1e3, 20*log10(max(syn(iz,:),1e-4)), 'LineWidth', 1.4); hold on;
plot(M.x*1e3, 20*log10(max(foc(iz,:),1e-4)), '--', 'LineWidth', 1.2);
plot(M.x*1e3, 20*log10(max(pw1(iz,:),1e-4)), ':', 'LineWidth', 1.2);
xlabel('x (mm)'); ylabel('dB (own max)'); grid on;
legend(sprintf('%d-angle sum', P.nAngles), 'focused', 'single PW');
title('lateral at z = 30 mm');

% (c) 纵向切片：合成场在焦点处出现峰值，平面波没有
subplot(2,2,4);
plot(M.z*1e3, 20*log10(max(syn(:,ix),1e-4)), 'LineWidth', 1.4); hold on;
plot(M.z*1e3, 20*log10(max(foc(:,ix),1e-4)), '--', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('dB (own max)'); grid on;
legend(sprintf('%d-angle sum', P.nAngles), 'focused'); title('axial at x = 0');
save_fig(fig, fullfile(F, 'F5_field_synthetic.png')); close(fig);

% --------------------------------------------------------------------- F6
% 一次点目标发射的原始通道数据。跨阵元的双曲线到达轨迹正是 DAS 要抵消的延时，
% 也是“为什么延时模型必须对”最直观的一张图。
S = load(fullfile(P.dataDir, 'pw', sprintf('P3_a%03d.mat', P.iZero)));
fig = figure('Position', [100 100 1000 420]);

% (a) 128 个通道画成图：横轴时间，纵轴阵元
subplot(1,3,[1 2]);
tt = (0:size(S.ch,1)-1)/P.fs*1e6;                 % 记录时间，µs
imagesc(tt, 1:P.N, 20*log10(abs(double(S.ch)).'/max(abs(double(S.ch(:))))));
axis xy; colormap(gca, gray); caxis([-60 0]);
xlabel('t (\mus)'); ylabel('element'); title('128-channel raw RF (\alpha = 0^\circ)');

% (b) 三条单通道波形：中心阵元最先到达
subplot(1,3,3);
for e = [16 64 112]
    plot(tt, double(S.ch(:,e))/max(abs(double(S.ch(:)))), 'LineWidth', 0.8); hold on;
end
xlim([30 50]); grid on; xlabel('t (\mus)'); ylabel('normalised');
legend('el 16','el 64','el 112'); title('three channels');
save_fig(fig, fullfile(F, 'F6_channel_rf.png')); close(fig);

% --------------------------------------------------------------------- F7
% 取包络之前的射频：PW 与 C11。两者并排，相干复合的效果才看得出来：单角度的
% 斜向条纹被压掉、散斑变细。这也是“求和发生在射频域而不是包络上”的直观证据。
fig = figure('Position', [100 100 820 380]);
ks = [1 2];                                       % PW、C11
for q = 1:2
    subplot(1,2,q);
    rr = double(R.P1.rf(:,:,ks(q)));
    rr = rr / max(abs(rr(:)));
    imagesc(P.x*1e3, P.z*1e3, 20*log10(max(abs(rr),1e-3)));
    axis image; colormap(gca, gray); caxis([-40 0]);
    xlim([-15 15]); ylim([18 45]);                % 放大焦点附近
    xlabel('x (mm)'); ylabel('z (mm)');
    title(sprintf('%s (RF before envelope)', R.names{ks(q)}));
end
save_fig(fig, fullfile(F, 'F7_rf_images.png')); close(fig);

% --------------------------------------------------------------------- F8
% 三幅并排：PW、C11、常规聚焦。全场统一 60 dB 色标，囊肿对比度可以直接用眼睛
% 比：C11 的囊肿更黑、边界更清楚，而聚焦图像只在 30 mm 焦点附近才好。
nPanel = numel(R.names);
fig = figure('Position', [100 100 1080 470]);
for q = 1:nPanel
    subplot(1,3,q);
    imagesc(P.x*1e3, P.z*1e3, R.P2.db(:,:,q)); axis image;
    colormap(gca, gray); caxis([-P.drDB 0]);
    title(sprintf('%s (n=%d)', R.names{q}, R.nAngles(q)), 'FontSize', 9);
    xlabel('x (mm)'); ylabel('z (mm)');
end
subplot(1,3,3);                                   % 第三格：聚焦参照
imagesc(R.FOC.x*1e3, P.z*1e3, R.FOC.db2); axis image; colormap(gca, gray);
caxis([-P.drDB 0]); title('focused (128 lines)', 'FontSize', 9);
xlabel('x (mm)'); ylabel('z (mm)');
save_fig(fig, fullfile(F, 'F8_bmode_cysts.png')); close(fig);

% --------------------------------------------------------------------- F9
% 均匀散斑体模用同样的并排布局：C11 的散斑明显更细、更均匀，而聚焦图像
% 只在焦点附近才均匀。
fig = figure('Position', [100 100 1080 470]);
for q = 1:nPanel
    subplot(1,3,q);
    imagesc(P.x*1e3, P.z*1e3, R.P1.db(:,:,q)); axis image;
    colormap(gca, gray); caxis([-P.drDB 0]);
    title(sprintf('%s (n=%d)', R.names{q}, R.nAngles(q)), 'FontSize', 9);
    xlabel('x (mm)'); ylabel('z (mm)');
end
subplot(1,3,3);                                   % 第三格：聚焦参照
imagesc(R.FOC.x*1e3, P.z*1e3, R.FOC.db1); axis image; colormap(gca, gray);
caxis([-P.drDB 0]); title('focused (128 lines)', 'FontSize', 9);
xlabel('x (mm)'); ylabel('z (mm)');
save_fig(fig, fullfile(F, 'F9_bmode_uniform.png')); close(fig);

% --------------------------------------------------------------------- F10
% 分辨率随深度。聚焦曲线是“只在焦点附近好”的参照（45–55 mm 退化到 1.6 mm）；
% 复合曲线是平的——这正是平面波复合的意义所在。
fig = figure('Position', [100 100 900 380]);

% (a) 横向：点目标的 -10 dB 宽度
subplot(1,2,1);
plot(R.zPts*1e3, R.w10.'*1e3, 'o-', 'LineWidth', 1.2); hold on;
plot(R.zPts*1e3, R.w10F*1e3, 'k--s', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('-10 dB lateral width (mm)'); grid on;
legend([R.names, {'focused'}], 'Location', 'northwest', 'FontSize', 7);
title('Lateral resolution');

% (b) 轴向：本来就应该与方法无关（复合不改变脉冲带宽），曲线确实都是平的
subplot(1,2,2);
plot(R.zPts*1e3, R.w6.'*1e3, 'o-', 'LineWidth', 1.2); hold on;
plot(R.zPts*1e3, R.w6F*1e3, 'k--s', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('-6 dB axial width (mm)'); grid on;
title('Axial resolution');
save_fig(fig, fullfile(F, 'F10_resolution.png')); close(fig);

% --------------------------------------------------------------------- F11
% 三种方法的对比度与 CNR 用柱状图并列：PW、C11、常规聚焦。误差棒是 9 个囊肿的
% 标准差。本实验只有一种复合设置，所以不再画“随角度数变化”的曲线。
fig = figure('Position', [100 100 900 380]);
names3 = [R.names, {'focused'}];
xb = (1:numel(names3)).';
Cc  = [arrayfun(@(s) s.Cmean, R.contrast),  R.contrastF.Cmean].';
Ce  = [arrayfun(@(s) s.Cstd,  R.contrast),  R.contrastF.Cstd].';
Cn  = [arrayfun(@(s) s.Cmean, R.contrastN), R.contrastFN.Cmean].';
CnE = [arrayfun(@(s) s.Cstd,  R.contrastN), R.contrastFN.Cstd].';

% (a) 对比度：无噪声与含噪声并列
subplot(1,2,1);
bar(xb, [Cc Cn], 'grouped'); hold on;
errorbar(xb-0.15, Cc, Ce, 'k.');  errorbar(xb+0.15, Cn, CnE, 'k.');
set(gca, 'XTick', xb, 'XTickLabel', names3);
ylabel('|C| (dB)'); grid on; legend('noise-free','with noise','Location','northwest');
title('Cyst contrast');

% (b) CNR：同样三组
subplot(1,2,2);
Cr  = [arrayfun(@(s) s.CNRmean, R.contrast), R.contrastF.CNRmean].';
CrE = [arrayfun(@(s) s.CNRstd,  R.contrast), R.contrastF.CNRstd].';
bar(xb, Cr); hold on; errorbar(xb, Cr, CrE, 'k.');
set(gca, 'XTick', xb, 'XTickLabel', names3);
ylabel('CNR'); grid on; title('CNR');
save_fig(fig, fullfile(F, 'F11_contrast.png')); close(fig);

% --------------------------------------------------------------------- F12
% 相对 SNR 随深度。灰色点线是各档位理论的 10log10(n)；画出来一眼就能看出实测
% 差了多少（差在哪、为什么，见报告 6.5）。
fig = figure('Position', [100 100 700 420]);
plot(R.zdep*1e3, R.snrRel.', 'o-', 'LineWidth', 1.2); hold on;
plot(R.zdep*1e3, R.snrRelF, 'k--s', 'LineWidth', 1.2);
for q = 1:numel(R.nAngles)
    plot(R.zdep([1 end])*1e3, 10*log10(R.nAngles(q))*[1 1], ':', 'Color', [0.6 0.6 0.6]);
end
xlabel('z (mm)'); ylabel('SNR relative to PW at 30 mm (dB)'); grid on;
legend([R.names, {'focused'}], 'Location', 'northeast', 'FontSize', 7);
title('SNR vs depth (relative)');
save_fig(fig, fullfile(F, 'F12_snr.png')); close(fig);

% --------------------------------------------------------------------- F13
% 散斑统计。实测强度分布在零附近比理想指数分布更尖，说明散斑没有完全发展
% （体模每个分辨单元只有几个散射子）。最右一幅改成三种方法的包络 μ/σ 柱状图，
% 并与 Rayleigh 参考值 1.91 对比。
sp = R.speck(1);
fig = figure('Position', [100 100 1000 380]);

% (a) 强度直方图，与 exp(-x) 对比（都是概率密度，所以直方图要除以区间宽度）
subplot(1,3,1);
bw = sp.histEdges(2) - sp.histEdges(1);
bar(sp.histEdges(1:end-1)+bw/2, sp.hist/bw, 1, 'FaceColor', [0.4 0.4 0.4]); hold on;
xx = linspace(0.01, 4, 200);
plot(xx, exp(-xx), 'r-', 'LineWidth', 1.2);       % 理想指数分布密度
xlabel('I / <I>'); ylabel('pdf'); grid on; legend('measured','exponential (ideal)');
title(sprintf('Intensity histogram, \\mu/\\sigma = %.2f (ideal 1.00)', sp.ratio));

% (b) 所有方法画在一起：形状随 n 变化很小
subplot(1,3,2);
for q = 1:numel(R.names)
    plot(0:numel(R.speck(q).hist)-1, R.speck(q).hist, '-', 'LineWidth', 1); hold on;
end
xlabel('bin'); ylabel('fraction per bin'); grid on; title('histograms, all methods');

% (c) 三种方法的包络 μ/σ，附 Rayleigh 参考线
subplot(1,3,3);
rA = [R.speck.ratioA, R.speckF.ratioA];
bar(1:3, rA); hold on;
plot([0.4 3.6], [1.91 1.91], 'k--');               % Rayleigh 参考
set(gca, 'XTick', 1:3, 'XTickLabel', [R.names, {'focused'}]);
ylabel('envelope \mu/\sigma'); grid on; title('speckle statistics');
save_fig(fig, fullfile(F, 'F13_speckle.png')); close(fig);

fprintf('figures written to %s\n', F);
end
