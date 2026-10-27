function make_figures()
%MAKE_FIGURES  由已保存的标定、声场与分析结果生成图 F1–F12。
%   在 run_tail / run_refresh 结尾调用。
%
%   输入（都在磁盘上，这里不重算任何仿真）
%     data/field/maps.mat     ：各角度的发射声场               （stage_fields）
%     data/res/results.mat    ：图像与指标                     （stage_analysis）
%     data/res/psfmap.mat     ：逐点散射子的 PSF 图             （run_psfmap）
%
%   例外：图 F2 画的是理想发射脉冲，直接由 params + pulse_waveform 生成，不读
%   任何数据文件——它展示的是“发射了什么”，不该随采集结果变化。
%
%   输出：figures/F1...F12.png，200 dpi。每张图保存后立即关闭，脚本可以无人值守
%   地跑完。
%
%   约定：所有坐标轴单位都是 mm；所有灰度 B 模式共用同一个 60 dB 色标
%   （P.displayRangeDB），这样各子图可以直接用眼睛比较。
%
%   图片清单
%   --------
%   F1  几何示意（阵列 + 三个波前）              —— 纯绘制
%   F2  理想发射脉冲的波形与频谱                 —— params（不经 Field II）
%   F3  单个 0° 平面波的发射声场                 —— maps
%   F4  九个偏转角的声场，显示波前倾斜           —— maps
%   F5  成像区域内的点扩散函数图（幅值/分辨率是否均匀） —— psfmap
%   F6  点目标的 128 通道原始射频                —— P3_a006（0°）
%   F7  取包络前的射频：PW / C11                  —— results
%   F8  囊肿体模 B 模式，各档位 + 聚焦           —— results
%   F9  散斑体模 B 模式，各档位 + 聚焦           —— results
%   F10 横向分辨率随深度                          —— results
%   F11 对比度与 CNR 随角度数                    —— results
%   F12 相对 SNR 随深度                          —— results

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);
M = load(fullfile(P.dataDir, 'field', 'maps.mat'));
R = load(fullfile(P.dataDir, 'res', 'results.mat'));

F = P.figureDir;
if ~exist(F, 'dir'), mkdir(F); end

% --------------------------------------------------------------------- F1
% 示意图，不是仿真结果：说明“偏转角 alpha 的平面波”在几何上是什么样子，
% 这样 F3/F4 的声场图才看得懂。
fig = figure('Position', [100 100 900 480]); hold on;

% 阵列：每四个阵元画一个小矩形（128 个全画就是一条黑杠）
for k = 1:4:P.nElements
    rectangle('Position', [(P.elementX(k)-P.elemWidth/2)*1e3, -0.35, P.elemWidth*1e3, 0.35], ...
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
% 理想发射脉冲（单程）的时域与频域。直接由 pulse_waveform 生成，不经过 Field II、
% 也不含接收端：这张图回答的是“发射的是什么”，不是“接收端量到什么”。
% 两程回波的实测值（约 4.1 MHz、半功率 62 %）在 T-A 表里另外给，两者不要混。
[imp, tp] = pulse_waveform(P);
ps = pulse_stats(imp, P.sampleRate);      % 与 E0 同一个测量函数，口径一致
fig = figure('Position', [100 100 900 340]);

% (a) 波形：高斯包络正弦，峰位在第 3 个载波周期；红虚线是包络
subplot(1,2,1);
env = abs(hilbert(imp)) / max(abs(imp));
plot(tp*1e6, imp/max(abs(imp)), 'LineWidth', 1.2); hold on;
plot(tp*1e6, env, 'r--', tp*1e6, -env, 'r--', 'LineWidth', 0.8);
xlabel('t (\mus)'); ylabel('normalised'); grid on; xlim([0 tp(end)*1e6]);
title(sprintf('Transmit pulse, f_c = %.2f MHz, \\sigma = %.2f T_c', ...
    P.pulseCarrierFreq/1e6, P.pulseSigmaCycles));

% (b) 频谱，同时标出 -6 dB 与半功率（-3 dB）两个电平的边沿
subplot(1,2,2);
plot(ps.f/1e6, ps.S, 'LineWidth', 1.2); hold on;
h6  = plot([ps.f1 ps.f2]/1e6,      [-6 -6], 'r-', 'LineWidth', 1.5);   % -6 dB
hhp = plot([ps.f1hp ps.f2hp]/1e6,  [-3 -3], 'm-', 'LineWidth', 1.5);   % 半功率
plot([ps.fc ps.fc]/1e6, [-45 0], 'k:');                                % 中心频率
xlim([1 9]); ylim([-45 2]); grid on;
xlabel('f (MHz)'); ylabel('dB');
% 图例必须显式绑定句柄：否则 legend 会依次取到频谱曲线和第一根横杠，颜色对不上
legend([h6 hhp], {'-6 dB (0.501)', 'half power (-3 dB, 0.708)'}, ...
    'Location', 'southwest');
title(sprintf('Transmit spectrum: -6 dB %.1f %% (half power %.1f %%)', ...
    100*ps.bw, 100*ps.bwHp));
save_fig(fig, fullfile(F, 'F2_pulse.png')); close(fig);

% --------------------------------------------------------------------- F3
% 0° 单平面波：没有焦点，近场最亮、随深度衰减（约 1/sqrt(z)）。这是与 F5 的
% 合成声场作对比的“之前”图。
i0 = find(P.angleIndex == 0);
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
% 成像区域到底均不均匀？在每个网格点上放一个散射子，走与正式采集完全相同的链
% 路（11 个角度发射 → 128 通道 → DAS → 射频域相加 → 包络），量该点的峰值幅度
% 与 -6 dB 宽度。
%
% 这张图取代了早先的“合成发射声场”：把 11 个角度的声场平移后相加，得到的是
% 一条窄条，那是“只补了轴上延时”留下的干涉残差，成像里并不存在这个加权。
% 真正被相加的是各角度的**图像**（DAS 已经逐像素对齐），所以成像区域是平的。
S = load(fullfile(P.dataDir, 'res', 'psfmap.mat'));
sx = S.xg*1e3;  sz = S.zg*1e3;                 % mm
[~, ix0] = min(abs(S.xg));                     % x = 0 那一列

% 深度方向的缓降是几何扩散（仿真里没有 TGC），不是加权缺陷。把它按深度取中值
% 当成“趋势”扣掉，剩下的残差才是“同一深度内各点是否均匀”的判据。
trend  = median(S.PdB, 2);                     % 每个深度的中值（列向量）
resid  = S.PdB - trend;                        % 残差图
aRes   = max(abs(resid(:)));

fig = figure('Position', [60 60 1180 780]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) 峰值幅度（绝对值）：整块区域只有一个很浅的上下渐变
nexttile;
imagesc(sx, sz, S.PdB); axis xy image; colormap(gca, parula);
cb = colorbar; cb.Label.String = 'dB re maximum';
hold on;
[XG, ZG] = meshgrid(sx, sz);
plot(XG(:), ZG(:), 'k+', 'MarkerSize', 3);     % 采样点（2 mm x 2.5 mm）
xlabel('x (mm)'); ylabel('z (mm)'); set(gca, 'FontSize', 9);
title(sprintf('(a) peak amplitude, %.1f .. %.1f dB', ...
    min(S.PdB(:)), max(S.PdB(:))), 'FontSize', 10);

% (b) 扣掉深度趋势之后的残差：这才是“各个点幅值是否均匀”的判据
nexttile;
imagesc(sx, sz, resid); axis xy image;
colormap(gca, div_cmap());                     % 蓝-白-红，0 在白色
cb = colorbar; cb.Label.String = 'dB re row median';
caxis(aRes*[-1 1]); xlabel('x (mm)'); ylabel('z (mm)'); set(gca, 'FontSize', 9);
title(sprintf('(b) after removing the depth trend: within %+.2f dB', aRes), ...
    'FontSize', 10);

% (c) 横向剖面：每个深度各自对齐到 x = 0，四条曲线几乎重合成一条直线
nexttile; hold on; grid on; box on;
zsel = [20 30 40 50];
for q = 1:numel(zsel)
    [~, iz] = min(abs(sz - zsel(q)));
    plot(sx, S.PdB(iz,:) - S.PdB(iz,ix0), 'o-', 'LineWidth', 1.2, 'MarkerSize', 3);
end
xlabel('x (mm)'); ylabel('dB re value at x = 0');
legend(arrayfun(@(v) sprintf('z = %d mm', v), zsel, 'UniformOutput', false), ...
       'Location', 'southwest', 'FontSize', 8);
title(sprintf('(c) lateral profiles, spread within a row %.2f dB', ...
    max(max(S.PdB,[],2) - min(S.PdB,[],2))), 'FontSize', 10);

% (d) 轴向：绝对幅度随深度缓降（几何扩散），中值线就是 (b) 扣掉的那条趋势
nexttile; hold on; grid on; box on;
plot(sz, S.PdB(:,ix0), 'o-', 'LineWidth', 1.3);
plot(sz, trend, 'k--', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('dB re maximum');
legend({'x = 0', 'row median'}, 'Location', 'southwest', 'FontSize', 8);
title(sprintf('(d) axial profile, %.1f dB over %d-%d mm', ...
    S.PdB(1,ix0) - S.PdB(end,ix0), round(sz(1)), round(sz(end))), 'FontSize', 10);
save_fig(fig, fullfile(F, 'F5_psf_map.png')); close(fig);

% --------------------------------------------------------------------- F6
% 一次点目标发射的原始通道数据。跨阵元的双曲线到达轨迹正是 DAS 要抵消的延时，
% 也是“为什么延时模型必须对”最直观的一张图。
S = load(fullfile(P.dataDir, 'pw', sprintf('P3_a%03d.mat', P.zeroAngleIndex)));
fig = figure('Position', [100 100 1000 420]);

% (a) 128 个通道画成图：横轴时间，纵轴阵元
subplot(1,3,[1 2]);
tt = (0:size(S.ch,1)-1)/P.sampleRate*1e6;                 % 记录时间，µs
imagesc(tt, 1:P.nElements, 20*log10(abs(double(S.ch)).'/max(abs(double(S.ch(:))))));
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
    imagesc(P.imageX*1e3, P.imageZ*1e3, 20*log10(max(abs(rr),1e-3)));
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
    imagesc(P.imageX*1e3, P.imageZ*1e3, R.P2.db(:,:,q)); axis image;
    colormap(gca, gray); caxis([-P.displayRangeDB 0]);
    title(sprintf('%s (n=%d)', R.names{q}, R.nAngles(q)), 'FontSize', 9);
    xlabel('x (mm)'); ylabel('z (mm)');
end
subplot(1,3,3);                                   % 第三格：聚焦参照
imagesc(R.FOC.x*1e3, P.imageZ*1e3, R.FOC.db2); axis image; colormap(gca, gray);
caxis([-P.displayRangeDB 0]); title('focused (128 lines)', 'FontSize', 9);
xlabel('x (mm)'); ylabel('z (mm)');
save_fig(fig, fullfile(F, 'F8_bmode_cysts.png')); close(fig);

% --------------------------------------------------------------------- F9
% 均匀散斑体模用同样的并排布局：C11 的散斑明显更细、更均匀，而聚焦图像
% 只在焦点附近才均匀。
fig = figure('Position', [100 100 1080 470]);
for q = 1:nPanel
    subplot(1,3,q);
    imagesc(P.imageX*1e3, P.imageZ*1e3, R.P1.db(:,:,q)); axis image;
    colormap(gca, gray); caxis([-P.displayRangeDB 0]);
    title(sprintf('%s (n=%d)', R.names{q}, R.nAngles(q)), 'FontSize', 9);
    xlabel('x (mm)'); ylabel('z (mm)');
end
subplot(1,3,3);                                   % 第三格：聚焦参照
imagesc(R.FOC.x*1e3, P.imageZ*1e3, R.FOC.db1); axis image; colormap(gca, gray);
caxis([-P.displayRangeDB 0]); title('focused (128 lines)', 'FontSize', 9);
xlabel('x (mm)'); ylabel('z (mm)');
save_fig(fig, fullfile(F, 'F9_bmode_uniform.png')); close(fig);

% --------------------------------------------------------------------- F10
% 横向分辨率随深度。聚焦只在 30 mm 的焦点处把宽度压到 0.77 mm，离开焦点就回到
% 0.95–1.04 mm（修正接收孔径后已不再“退化到 1.6 mm”，见
% docs/实验报告_讨论与附录.md 的 B12）；复合曲线是平的——这正是平面波复合的意义。
% 轴向分辨率不再画：它由脉冲长度 c*tau/2 决定，三种方法实测都是 0.19 mm
% （0.176–0.199 mm，起伏在读数不确定度以内），比不出差别（同文件 B14）。
fig = figure('Position', [100 100 620 430]);
plot(R.zPts*1e3, R.w10.'*1e3, 'o-', 'LineWidth', 1.2); hold on;
plot(R.zPts*1e3, R.w10F*1e3, 'k--s', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('-10 dB lateral width (mm)'); grid on;
legend([R.names, {'focused'}], 'Location', 'southwest', 'FontSize', 7);
title('Lateral resolution');
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
% 差了多少（差在哪、为什么，见报告 5.7）。
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

fprintf('figures written to %s\n', F);
end

% =========================================================================
function cm = div_cmap()
%DIV_CMAP  蓝-白-红发散色标：0 落在白色，正负一眼可分（残差图专用）。
%   MATLAB 没有内置发散色标，这里用两段线性渐变拼一个，只有 128 级，够用。
up = [linspace(0.15, 1, 64).' linspace(0.35, 1, 64).' linspace(0.75, 1, 64).'];  % 蓝 -> 白
dn = [linspace(1, 0.70, 64).' linspace(1, 0.10, 64).' linspace(1, 0.10, 64).'];  % 白 -> 红
cm = [up; dn];
end
