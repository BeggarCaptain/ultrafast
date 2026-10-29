function make_figures()
%MAKE_FIGURES  由已保存的标定、声场与分析结果生成图 F1–F12。
%   在 run_tail / run_refresh 结尾调用。
%
%   输入（都在磁盘上，这里不重算任何仿真）
%     data/field/maps.mat     ：各角度的发射声场               （stage_fields）
%     data/res/results.mat    ：图像与指标                     （stage_analysis）
%     data/res/psfmap.mat     ：逐点散射子的 PSF 图             （run_psfmap）
%     data/res/axisres.mat    ：轴上逐点单散射子实验             （run_axisres）
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
%   约定：换能器一律画在图片上方，z 由上至下增大
%   ------------------------------------------
%   这是超声图像的标准摆法（探头在上、深度向下），也便于把仿真图直接与临床/
%   文献图对照。实现上就是保持 imagesc 的默认 YDir = 'reverse'：不要再写
%   axis xy（它会把 y 轴翻正，变成深度向上）；需要在 plot 里画示意图时显式
%   set(gca,'YDir','reverse')。折线类子图（F2、F5(c)(d)、F10–F12）不是图像，
%   按各自的自变量画，不受这条约定约束。
%
%   图片清单
%   --------
%   F1  几何示意（阵列 + 三个波前）              —— 纯绘制
%   F2  理想发射脉冲的波形与频谱                 —— params（不经 Field II）
%   F3  单个 0° 平面波的发射声场                 —— maps
%   F4  六个偏转角的声场，显示波前倾斜           —— maps
%   F5  成像区域内的点扩散函数图（幅值/分辨率是否均匀） —— psfmap
%   F6  点目标：128 通道射频幅度 → 单角度 DAS → 11 角度复合 —— P3_a006（0°）
%   F7  取包络前的射频：PW / C11                  —— results
%   F8  囊肿体模 B 模式，各档位 + 聚焦           —— results
%   F9  散斑体模 B 模式，各档位 + 聚焦           —— results
%   F10 横向分辨率随深度（轴上 26 个深度逐点重扫） —— axisres
%   F11 囊肿对比度（三种方法并列）               —— results
%   F12 SNR 随深度（绝对刻度）                    —— results
%
%   为什么按这个顺序排（每张图回答一个问题）
%   --------------------------------------
%   报告的证据链是"发射 → 声场 → 图像 → 指标"，图也照这条链走：
%     F1/F2         发射的是什么（几何、脉冲）——后面所有结果的输入；
%     F3/F4/F14     声场长什么样、11 束覆盖到哪儿——解释分辨率与 SNR 为什么随
%                   位置变化，也是"复合只在公共照明区成立"的依据；
%     F5            成像区域内每个点的 PSF 是否一致——证明复合不是"只在焦点好"；
%     F6/F7         波束形成前后的原始信号——让读者相信图像确实由通道数据算出来，
%                   而不是某个中间量；
%     F8/F9         三种方法的 B 模式对照——直观印象；
%     F10–F12       三个指标的定量曲线——结论就是这三张。
%   这样排，读图的人不必先接受指标，可以顺着链自己看出结论。
%
%   为什么画图阶段只读磁盘、绝不重算
%   ------------------------------
%   如果画图时顺手重算声场或指标，一旦图与表对不上，就分不清是数据错了还是画错
%   了。这里全部只读：图是数据的确定性函数，数据是报告的确定性函数，出错时能
%   二分定位。
%
%   为什么统一色标与单位
%   ------------------
%   所有 B 模式图共用 60 dB 动态范围、所有坐标轴用 mm，子图之间才能直接用眼睛
%   比较；否则"看起来更亮"可能只是色标不同。色标的具体范围在各子图的 caxis
%   处注明。

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

text(-24, 4, '128 阵元线阵，间距 0.33 mm', 'FontSize', 9);
text(6, 55, '\alpha = +10\circ', 'FontSize', 9);
text(6, 57, '\alpha = 0\circ', 'FontSize', 9);
text(-24, 55, '\alpha = -10\circ', 'FontSize', 9);
xlabel('x (mm)'); ylabel('z (mm)'); axis([-25 25 -2 62]); grid on;
set(gca, 'YDir', 'reverse');      % 换能器在上：深度从上往下增大（全项目统一）
title('平面波发射几何');
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
xlabel('t (\mus)'); ylabel('归一化'); grid on; xlim([0 tp(end)*1e6]);
title(sprintf('发射脉冲：f_c = %.2f MHz，\\sigma = %.2f T_c', ...
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
legend([h6 hhp], {'−6 dB（幅度 0.501）', '半功率（−3 dB，0.708）'}, ...
    'Location', 'southwest');
title(sprintf('发射频谱：−6 dB %.1f %%，半功率 %.1f %%', ...
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
title('发射声压，\alpha = 0\circ（dB）');

% (b) 30 mm 处横向切片：基本是平的，即没有横向约束
subplot(2,2,2);
[~, iz] = min(abs(M.z - 30e-3));
plot(M.x*1e3, 20*log10(max(mp(iz,:), 1e-4)), 'LineWidth', 1.2);
xlabel('x (mm)'); ylabel('dB'); grid on; title('z = 30 mm 处横向剖面');

% (c) 轴线上的纵向切片：单调衰减，没有焦峰
subplot(2,2,4);
[~, ix] = min(abs(M.x));
plot(M.z*1e3, 20*log10(max(mp(:,ix), 1e-4)), 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('dB'); grid on; title('x = 0 处轴向剖面');
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
%
% 读图要点
%   2 x 2：上排是峰值幅度，下排是横向 −6 dB 宽度；每排都是“地图 + 扣掉深度趋势的
%   残差”。(a) 整块只有上下很浅的渐变（几何扩散，仿真里没有 TGC）；(b) 扣掉深度
%   趋势后残差在 ±1 dB 量级，且不出现随深度收敛的条纹；(c) 横向宽度 0.71–0.78 mm；
%   (d) 残差在 ±0.05 mm 以内，说明同一深度内横向几乎一样宽。
%   如果幅度残差 (b) 里出现一条沿深度收窄的亮带，那就是“复合合成了一个发射焦点”的
%   假象——正是 B13 指出的那种错误加权。现在没有，说明复合只是把各角度的图像叠加，
%   没有制造新的焦区。
%
%   颜色约定：**暖色 = 指标更好**。幅度越大越好，所以 (a)(b) 用常规 parula/发散色标
%   （亮 = 强 = 好）；横向分辨率越小越好，所以 (c)(d) 的色标整个翻过来（窄 = 亮 =
%   好，宽 = 暗 = 差）。否则上下两排同一个颜色在说相反的话，看的人得来回换算。
%   排成 2 x 2 是为了让同类的“地图 + 残差”左右相邻、幅度与分辨率上下相邻；成像区域
%   是 20 x 35 mm 的竖长条，四幅按 1:1 比例正好占满一个竖幅（缩到报告页宽后每幅约
%   3.3 cm 宽）。
S = load(fullfile(P.dataDir, 'res', 'psfmap.mat'));
sx = S.xg*1e3;  sz = S.zg*1e3;                 % mm
[~, ix0] = min(abs(S.xg));                     % x = 0 那一列

% 深度方向的缓降是几何扩散（仿真里没有 TGC），不是加权缺陷。把它按深度取中值
% 当成“趋势”扣掉，剩下的残差才是“同一深度内各点是否均匀”的判据。
trend  = median(S.PdB, 2);                     % 每个深度的中值（列向量）
resid  = S.PdB - trend;                        % 残差图
aRes   = max(abs(resid(:)));

% 图幅比例是为报告页留的：四幅 1:1 的竖长条按 2 x 2 排，内容本身约 40:70（宽:高），
% 直接把图幅设成这个比例的话，缩到 A4 页宽（14.8 cm）后有 22 cm 高、图注会被挤到
% 下一页；这里把图幅放宽到 1150:1500，多出来的宽度变成两列之间的留白，图高降到
% 19.3 cm，图注就能跟在同一页。
fig = figure('Position', [60 60 1150 1500]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) 峰值幅度（绝对值）：整块区域只有一个很浅的上下渐变
nexttile;
imagesc(sx, sz, S.PdB); axis image; set(gca, 'YDir', 'reverse');   % 换能器在上
colormap(gca, parula);
cb = colorbar; cb.Label.String = 'dB（相对最大值）';
hold on;
[XG, ZG] = meshgrid(sx, sz);
plot(XG(:), ZG(:), 'k+', 'MarkerSize', 3);     % 采样点（2 mm x 2.5 mm）
xlabel('x (mm)'); ylabel('z (mm)'); set(gca, 'FontSize', 9);
title(sprintf('(a) 峰值幅度：%.1f … %.1f dB', ...
    min(S.PdB(:)), max(S.PdB(:))), 'FontSize', 10);

% (b) 扣掉深度趋势之后的残差：这才是“各个点幅值是否均匀”的判据
nexttile;
imagesc(sx, sz, resid); axis image; set(gca, 'YDir', 'reverse');   % 换能器在上
colormap(gca, div_cmap());                     % 蓝-白-红，0 在白色
cb = colorbar; cb.Label.String = 'dB（相对该深度中值）';
caxis(aRes*[-1 1]); xlabel('x (mm)'); ylabel('z (mm)'); set(gca, 'FontSize', 9);
title(sprintf('(b) 扣除深度趋势后的残差：起伏在 %+.2f dB 以内', aRes), ...
    'FontSize', 10);

% (c) 横向 −6 dB 宽度（同一批点）：整块区域几乎不变，没有随深度收敛
WXmm   = S.WX * 1e3;                           % mm
wTrend = median(WXmm, 2);                      % 每个深度的中值
WXres  = WXmm - wTrend;                        % 残差
aW     = max(abs(WXres(:)));
nexttile;
imagesc(sx, sz, WXmm); axis image; set(gca, 'YDir', 'reverse');   % 换能器在上
colormap(gca, flipud(parula));                 % 翻过来：窄 = 亮 = 好
cb = colorbar; cb.Label.String = 'mm（越暗越宽）';
caxis([0.70 0.79]); hold on;
plot(XG(:), ZG(:), 'k+', 'MarkerSize', 3);     % 采样点，与 (a) 同一批
xlabel('x (mm)'); ylabel('z (mm)'); set(gca, 'FontSize', 9);
title(sprintf('(c) −6 dB 横向宽度：%.2f … %.2f mm', ...
    min(WXmm(:)), max(WXmm(:))), 'FontSize', 10);

% (d) 扣掉逐深度中值后的残差：同一深度内横向是否处处一样宽
nexttile;
imagesc(sx, sz, WXres); axis image; set(gca, 'YDir', 'reverse');  % 换能器在上
colormap(gca, flipud(div_cmap()));             % 与 (b) 同一套，但翻过来：变窄 = 红
cb = colorbar; cb.Label.String = 'mm（红 = 更窄）';
caxis(aW*[-1 1]); xlabel('x (mm)'); ylabel('z (mm)'); set(gca, 'FontSize', 9);
title(sprintf('(d) 扣除深度趋势后的残差：起伏在 ±%.3f mm 以内', aW), ...
    'FontSize', 10);
save_fig(fig, fullfile(F, 'F5_psf_map.png')); close(fig);

% --------------------------------------------------------------------- F6
% 三个台阶并排：通道数据 → 单角度成像 → 11 角度复合。
%   (a) 同一帧的 128 通道射频幅度，跨阵元的双曲线到达轨迹正是 DAS 要抵消的延时；
%   (b) 这组通道数据经 DAS + 包络后的 B 模式，双曲线被拉平、点目标聚成亮点；
%   (c) 再加 10 个角度做相干复合（C11）：单角度残留的斜向旁瓣被压掉，点更紧。
%
% 三幅共用同一套坐标与比例：横轴都是 mm，纵轴都是深度 mm。(a) 的横轴是阵元位置
% （±21.12 mm），纵轴由记录时刻换算成往返深度 z = c*t/2，所以“换能器在上、深度
% 向下”的读法与 (b)(c) 完全一致，同一个点目标在三幅里位于同一高度。
%
% 记录窗口从第一个可能回波起（t0 = 19.47 µs，即 z = 15 mm），到 60.9 mm 为止，
% 所以三幅都只画 z = 15–60 mm；体模最近的散射子就在 15 mm。
S = load(fullfile(P.dataDir, 'pw', sprintf('P3_a%03d.mat', P.zeroAngleIndex)));
ch = double(S.ch) / max(abs(double(S.ch(:))));          % 全局归一化到 0 dB
xe = P.elementX * 1e3;                                   % 阵元位置 [mm]
zch = (S.t0 + (0:size(ch,1)-1).'/P.sampleRate) * P.soundSpeed/2 * 1e3;   % 往返深度 [mm]
axX = [-21.2 21.2];  axZ = [15 60];                      % 三幅共用

fig = figure('Position', [60 60 1750 720]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) 128 通道原始射频：换能器在上，纵轴是往返深度
%
% 画的是射频的幅度（取包络后再显示）。载波 4.5 MHz 在 30 MHz 采样下只有 6.7 个
% 采样/周期，而 45 mm 的画面只有几百像素高——直接画 |RF| 会把载波抽样成一片网格
% 状的摩尔纹，双曲线的到达轨迹反而看不清。包络不丢任何延时信息（它正是 DAS 之后
% 要叠加的那个量），只是把每个通道的振荡去掉。
env = abs(hilbert(ch));                                  % 每通道的射频包络
nRow = 700;                                              % 子图大约的像素高度
edge = round(linspace(1, size(env,1) + 1, nRow + 1));
envD = zeros(nRow, size(env,2), 'single');
zD   = zeros(nRow, 1);
for r = 1:nRow
    ii = edge(r):edge(r+1)-1;
    envD(r,:) = single(mean(env(ii,:), 1));
    zD(r)     = mean(zch(ii));
end
nexttile;
imagesc(xe, zD, 20*log10(max(envD, 1e-4)));      % envD 是 [nRow x 128]：行=深度、列=阵元
axis image; set(gca, 'YDir', 'reverse');                 % 换能器在上
colormap(gca, gray); caxis([-60 0]);
xlabel('阵元位置 x (mm)'); ylabel('z = ct/2 (mm)');
title('(a) 128 通道原始射频（\alpha = 0^\circ）', 'FontSize', 10);
xlim(axX); ylim(axZ);

% (b)(c) 同一组通道数据：单角度成像 与 11 角度复合。两幅都用 B 模式、60 dB，
% 各自按自身最大值归一化（to_db 的口径），所以看的是形状与旁瓣，不是绝对亮度。
ttl = {'(b) 0^\circ 单角度 DAS 成像（PW）', '(c) 11 角度复合（C11）'};
for k = 1:2
    nexttile;
    imagesc(P.imageX*1e3, P.imageZ*1e3, R.P3.db(:,:,k));
    axis image; set(gca, 'YDir', 'reverse');             % 换能器在上
    colormap(gca, gray); caxis([-P.displayRangeDB 0]);
    xlabel('x (mm)'); ylabel('z (mm)');
    title(ttl{k}, 'FontSize', 10);
    xlim(axX); ylim(axZ);
end
cb = colorbar; cb.Label.String = 'dB';                   % 三幅同一色标（各自 0 dB）
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
    title(sprintf('%s（取包络前的射频）', R.names{ks(q)}));
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
    title(sprintf('%s（%d 次发射）', R.names{q}, R.nAngles(q)), 'FontSize', 9);
    xlabel('x (mm)'); ylabel('z (mm)');
end
subplot(1,3,3);                                   % 第三格：聚焦参照
imagesc(R.FOC.x*1e3, P.imageZ*1e3, R.FOC.db2); axis image; colormap(gca, gray);
caxis([-P.displayRangeDB 0]); title('常规聚焦（128 线）', 'FontSize', 9);
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
    title(sprintf('%s（%d 次发射）', R.names{q}, R.nAngles(q)), 'FontSize', 9);
    xlabel('x (mm)'); ylabel('z (mm)');
end
subplot(1,3,3);                                   % 第三格：聚焦参照
imagesc(R.FOC.x*1e3, P.imageZ*1e3, R.FOC.db1); axis image; colormap(gca, gray);
caxis([-P.displayRangeDB 0]); title('常规聚焦（128 线）', 'FontSize', 9);
xlabel('x (mm)'); ylabel('z (mm)');
save_fig(fig, fullfile(F, 'F9_bmode_uniform.png')); close(fig);

% --------------------------------------------------------------------- F10
% 横向分辨率随深度：x = 0 轴上每 2 mm 放一个散射子，一次只放一个，逐个重扫
% （data/res/axisres.mat，26 个深度 10–60 mm）。旧版只有点目标体模的 6 个深度，
% 画不出曲线形状，现在改用这条逐点实验。
%
% 三条 2 mm 曲线，外加一条聚焦的 0.5 mm 细扫（浅灰）：
%   * PW   0.83–0.87 mm，随深度很平；
%   * C11  0.74–0.76 mm 的平线，比 PW 窄 11–13 %——这正是复合的意义；
%   * 聚焦 呈"V"形，最低 0.59–0.65 mm 落在 30–32 mm（发射焦点 30 mm），
%          焦外升到 0.90–0.95 mm。
%
% 细扫那条线不是插值，也不是"更精确的版本"：灰线上每个点同样是一次独立的
% 单散射子实验（每 0.5 mm 放一个、各自成像、各自量宽度），与三角标记在重合的
% 26 个深度上逐位相同。它只解决一件事——沿深度的网格够不够密：场本身在焦区
% 有约 1.8 mm 周期的轴向纹波（矩形发射孔径的边缘波，轴上峰值声压起伏 4–6 dB），
% 点目标落在波峰/波谷上宽度最多差 0.1 mm。2 mm 的 26 个点分辨不出这个尺度，
% 相邻两点连成直线就会把 25 mm 的 0.90 抹成"24→26 mm 的缓降"（0.80→0.72）；
% 表 T-B 的 25/35 mm 也正好落在两点之间，只有细扫才在图上找得到对应位置。
% 注意用词：这不是"混叠"（每个点都是精确的独立测量，没有重建），只是曲线在
% 两点之间未知。
%
% 读图要点：两条平线 + 一条带纹波的 V。若聚焦曲线在焦外出现 1.5 mm 以上的尖峰，
% 说明接收孔径又被发射窗截断了（B12 的历史问题复现）。
A = load(fullfile(P.dataDir, 'res', 'axisres.mat'));
fig = figure('Position', [100 100 680 470]);
hold on; box on;
% 颜色显式给出：先画的细扫线如果不指定颜色，会把默认色序的第一格占掉，
% 后面两条曲线就会错位成橙/黄。
cPW = [0 0.447 0.741];  cC11 = [0.850 0.325 0.098];
% 聚焦的 0.5 mm 细扫先画，压在标记下面
h0 = plot(A.zgFine*1e3, A.wFocFine*1e3, '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.7);
h1 = plot(A.zg*1e3, A.wPW*1e3,  'o-', 'Color', cPW,  'LineWidth', 1.3, 'MarkerSize', 4);
h2 = plot(A.zg*1e3, A.wC11*1e3, 's-', 'Color', cC11, 'LineWidth', 1.3, 'MarkerSize', 4);
h3 = plot(A.zg*1e3, A.wFoc*1e3, 'k--^', 'LineWidth', 1.1, 'MarkerSize', 5, ...
          'MarkerFaceColor', 'w');
xlim([10 60]); ylim([0.57 1.0]);
xlabel('z (mm)'); ylabel('−6 dB 横向宽度（mm）'); grid on;
% 图例放在坐标区外：焦区那段曲线（0.59–0.95）恰好占满坐标区左下角，
% 图例放里面会把最关键的谷底挡住。
legend([h1 h2 h3 h0], {'单平面波（PW）', '11 角度复合（C11）', ...
        '常规聚焦（2 mm 采样）', '常规聚焦（0.5 mm 细扫）'}, ...
       'Location', 'southoutside', 'Orientation', 'horizontal', ...
       'NumColumns', 2, 'FontSize', 7.5);
title('横向分辨率随深度（轴上逐点单散射子）');
save_fig(fig, fullfile(F, 'F10_resolution.png')); close(fig);

% --------------------------------------------------------------------- F11
% 三种方法的囊肿对比度用柱状图并列：PW、C11、常规聚焦，无噪声与含噪声成对。
% 误差棒是 9 个囊肿的标准差。本实验只有一种复合设置，所以不画“随角度数变化”的曲线。
%
% 读图要点：C11 的 |C| 比 PW 高约 10.5 dB（23.3 对 12.8）——复合压低了旁瓣
% 泄漏，囊肿内部更"黑"；聚焦的均值与 C11 接近（22.2）但标准差大 2.4 倍，因为
% 它把对比度堆在焦区。含噪柱子比无噪低 0.9–2.1 dB（噪声底把囊肿内部填起来
% 了），这正是两套对比度要分开报的原因。
%
% CNR 不出图：三组都在 0.58–0.63，差别来自散斑纹理而不是囊的可见度（方向还与
% |C| 相反），单独画一幅只会诱导读者把它读成质量排序（原因见附录 B23）。
% 数值仍列在表 T-C。
fig = figure('Position', [100 100 560 380]);
names3 = [R.names, {'单焦点'}];
xb = (1:numel(names3)).';
Cc  = [arrayfun(@(s) s.Cmean, R.contrast),  R.contrastF.Cmean].';
Ce  = [arrayfun(@(s) s.Cstd,  R.contrast),  R.contrastF.Cstd].';
Cn  = [arrayfun(@(s) s.Cmean, R.contrastN), R.contrastFN.Cmean].';
CnE = [arrayfun(@(s) s.Cstd,  R.contrastN), R.contrastFN.Cstd].';

bar(xb, [Cc Cn], 'grouped'); hold on;
errorbar(xb-0.15, Cc, Ce, 'k.');  errorbar(xb+0.15, Cn, CnE, 'k.');
set(gca, 'XTick', xb, 'XTickLabel', names3);
ylabel('|C|（dB）'); grid on; legend('无噪声','含噪声','Location','northwest');
title('囊肿对比度');
save_fig(fig, fullfile(F, 'F11_contrast.png')); close(fig);

% --------------------------------------------------------------------- F12
% 绝对 SNR 随深度。灰点线是"单平面波 30 mm 处"这条刻度锚点加上理论增益
% 10log10(n)，画出来一眼就能看出实测比理想相干复合低多少（差在哪、为什么，
% 见报告 5.7）。
%
% 读图要点：三条曲线在 30 mm 附近约为 20.4 / 29.2 / 31.3 dB；C11 比它自己那条
% 灰线（锚点 + 10.4 dB = 30.8 dB）低约 1.6 dB，且随深度基本平；单焦点在 30 mm
% 最高，向两侧快速下降——它的"高"来自把发射能量集中在焦点，不是复合增益，所以
% 形状完全不同。绝对刻度的含义：σ_e 按"PW 在 30 mm 处 20 dB"标定，所以这些数
% 不是临床可比的绝对灵敏度；有意义的量是方法之间的差（C11−PW = 8.8 dB）。
fig = figure('Position', [100 100 700 420]);
plot(R.zdep*1e3, R.snr.', 'o-', 'LineWidth', 1.2); hold on;
plot(R.zdep*1e3, R.snrF, 'k--s', 'LineWidth', 1.2);
for q = 1:numel(R.nAngles)
    yb = R.snrAbsPW30 + 10*log10(R.nAngles(q));      % 锚点 + 理论增益
    plot(R.zdep([1 end])*1e3, yb*[1 1], ':', 'Color', [0.6 0.6 0.6]);
end
xlabel('z (mm)'); ylabel('SNR（dB）'); grid on;
legend([R.names, {'单焦点'}], 'Location', 'northeast', 'FontSize', 7);
title('SNR 随深度');
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
