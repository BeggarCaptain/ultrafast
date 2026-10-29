function stage_fields(P, force)
%STAGE_FIELDS  各平面波的发射峰值声压图。
%
%   STAGE_FIELDS(P) 用 calc_hp 在 0.5 mm 网格（x = -20..20 mm、z = 5..60 mm，
%   共 81 x 111 点）上计算每个角度的发射峰值声压 -> peakAll(:,:,i)，存
%   data/field/maps.mat。
%
%   STAGE_FIELDS(P, true) 忽略缓存强制重算（改了声场网格、或想核对结果时用）。
%
%   这份数据有两个用途
%     * F3/F4：单角度与偏转角声场的形状（波前倾斜、没有发射焦点）；
%     * F14   ：11 个角度的幅度和 sum_alpha |p_alpha|，即 DAS 之后图像实际看到
%               的发射加权。
%
%   为什么要亚样本平均（修正见附录 B22）
%   -----------------------------------
%   calc_hp 返回的 p(t) 采样在 1/fs = 33.3 ns 的网格上，而 0.30 mm 阵元的空间
%   冲激响应边缘不到 1 us，相对这个网格就是一串台阶。某个格点的回波落在网格的
%   哪个位置，由它的到达时间决定；到达时间每加深 c/fs = 51.33 um 正好推移一个
%   样本。于是幅度上出现一条周期 51.33 um、峰峰 1-2 dB 的伪迹——它不是声场结构：
%   把同一个点的输出采样率提到 60/120/240 MHz，剖面立刻变单调，与 240 MHz 参考
%   的偏差从 1.07 dB 依次降到 0.36 / 0.13 / 0 dB。
%
%   51.33 um 比场图网格小两个数量级，直接采样必然混叠：1 mm 一格里含 19.48 个
%   周期，小数部分 0.48，于是相邻行反相——早期 1 mm 网格上图里的"横纹"就是这么
%   来的。现在每个格点在 N 个亚样本偏移 z + k*c/(N*fs)（k = 0..N-1）上各算一次
%   再取平均，正好覆盖一个完整周期（与偏转角无关：平移 z 就是把到达时间平移
%   1/(N*fs)），把相位依赖平均掉，只留真实声场。
%
%   为什么网格用 0.5 mm
%   -------------------
%   这张图要看的是"波前倾斜成什么样、11 束覆盖到哪"，尺度是厘米；0.5 mm 的格子
%   在 40 x 55 mm 的画面上已经细到看不出台阶（1 mm 时 flat shading 本身就是一条
%   一条 1 mm 高的色带），再加密一倍计算量翻两番，换来的只是更细的等值线边缘。
%   要量 PSF 宽度是另一件事，那由 run_psfmap 在每个点周围开 0.05 mm 的局部小图
%   来做（见 psf_map）。
%
%   代价：81 x 111 x 11 个点 x 5 个亚样本，约 80 s（1 mm 网格、单个亚样本时约 4 s）。
%
%   为什么只取 max(abs(hp))、不记峰值时刻
%   ------------------------------------
%   calc_hp 返回的是每个网格点上的声压**波形**。这里只需要"这个点被照得多亮"，
%   幅度就够；时间信息由 DAS 的延时模型负责，场图不参与成像计算。也正因如此，
%   记录长度与时间基准都不必对齐，省掉一堆簿记。
%
%   为什么每个角度单独算，而不是先合成再算
%   ------------------------------------
%   声场满足线性叠加，但"合成"这一步在成像里是逐像素补相位之后才做的（F14b、
%   psf_map）。若先把 11 个角度的场在空间上平移相加，得到的是只补了一半延时的
%   干涉残差，不是成像里存在的任何量——这条弯路已经走过，见附录 B13。所以这里
%   老老实实存 11 个角度各自的结果，合成留给读图的人按正确方式做。
%
%   为什么不再保存“11 个角度声场平移后相加”的合成场
%   ----------------------------------------------
%   那种合成场只把每个角度的延时补到轴上一个常数（z0*cos(alpha)/c），横向相位
%   项 x*sin(alpha) 原样留着，于是叠加结果只是一条窄条：宽度约
%   c*T/(2*sin(delta_alpha))、轴向尺度约 lambda/(1-cos(delta_alpha))。
%   它不是成像里存在的任何量——真实成像在 DAS 里逐像素补掉 (x*sin+z*cos)/c，
%   图像看到的发射加权是幅度相加 sum|p|（见 F14），那条窄条只是“只补了一半
%   延时”留下的干涉残差。要检验成像区域，正确做法是逐点成像（run_psfmap → F5）。
%
%   为什么用 calc_hp 而不是 calc_scat_multi
%   ---------------------------------------
%   calc_hp 直接给出指定点上的发射声压波形，用来画“场”；calc_scat_multi 返回的
%   是加了散射子、接收延时与脉冲回波之后的通道信号，反推场反而更绕。两者共用
%   同一套阵列与激励，所以 calc_hp 的结果就是这个系统真实的发射场。
%
%   缓存：maps.mat 里存一份输入指纹（src/result_stamp.m），参数、网格、亚样本数与
%   计算逻辑都没变就直接返回。这是确定性计算——同样的输入必然给出逐位相同的
%   输出，重算没有信息量。详见 src/tag_match.m。
%
%   输出用 single 精度；原始声压记录不保存（每次调用约 85 MB）。

fprintf('\n===== transmit field maps =====\n');
if nargin < 2, force = []; end
force = force_flag(force);
outdir = fullfile(P.dataDir, 'field');
if ~exist(outdir, 'dir'), mkdir(outdir); end

% 指纹要含网格与亚样本数：fieldMapX/Z/fieldMapSubSamples 不在 data_stamp 里
% （改它们不该触发重跑采集），但它们决定了 maps.mat 的内容。
tag = result_stamp(P, struct('x', P.fieldMapX, 'z', P.fieldMapZ, ...
                             'nsub', P.fieldMapSubSamples));
out = fullfile(outdir, 'maps.mat');
if ~force && tag_match(out, tag)
    fprintf('field maps: 参数未变，跳过（%s）\n', out);
    return;
end

x = P.fieldMapX;  z = P.fieldMapZ;                       % 0.5 mm 网格
[Xg, Zg] = meshgrid(x, z);
nSub = P.fieldMapSubSamples;                             % 亚样本偏移的个数
dzSub = P.soundSpeed / P.sampleRate / nSub;              % 每个偏移 [m]，5 个正好一个 c/fs

Th = f2_setup(P);
peakAll = zeros(numel(z), numel(x), P.nAngles, 'single');
t0 = tic;

for i = 1:P.nAngles
    a = P.angles(i);
    d = P.elementX * sin(a) / P.soundSpeed;                 % 与采集时完全相同的发射延时
    xdc_center_focus(Th, [0 0 0]);
    xdc_focus_times(Th, 0, d);
    xdc_apodization(Th, 0, ones(1, P.nElements));    % 时间线 0 s 起生效（不是通道编号）

    % N 个亚样本偏移累加：每个偏移只改 z，到达时间就平移 1/(N*fs)，
    % 5 个点覆盖一个完整的 c/fs 周期，伪迹的相位依赖被平均掉。
    acc = zeros(numel(Zg), 1);
    for k = 0:nSub-1
        pts = [Xg(:), zeros(numel(Xg), 1), Zg(:) + k*dzSub];   % calc_hp 要 N x 3
        hp  = calc_hp(Th, pts);                      % [Nt x npts] 声压
        % 每个点只关心“峰值有多大”，不关心它出现在哪个时刻：记录长度与时间基准
        % 都不需要对齐，取 max(abs()) 即可。
        acc = acc + max(abs(hp), [], 1).';
    end
    peakAll(:,:,i) = single(reshape(acc / nSub, numel(z), numel(x)));
end

S = struct('x', x, 'z', z, 'peakAll', peakAll, 'nsub', nSub, 'tag', tag);
save(out, '-struct', 'S', '-v7.3');   % 逐角度声场，用 v7.3
fprintf('saved field maps (%d angles, %d x %d grid, %d sub-samples, %.1f s)\n', ...
        P.nAngles, numel(x), numel(z), nSub, toc(t0));
end
