function stage_fields(P)
%STAGE_FIELDS  各平面波的发射峰值声压图。
%
%   STAGE_FIELDS(P) 用 calc_hp 在 1 mm 网格（x = -20..20 mm、z = 5..60 mm，共
%   41 x 56 点）上计算每个角度的发射峰值声压 -> peakAll(:,:,i)，存
%   data/field/maps.mat。
%
%   这份数据有两个用途
%     * F3/F4：单角度与偏转角声场的形状（波前倾斜、没有发射焦点）；
%     * F14   ：11 个角度的幅度和 sum_alpha |p_alpha|，即 DAS 之后图像实际看到
%               的发射加权。
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
%   输出用 single 精度；原始声压记录不保存（每个角度约 44 MB）。

fprintf('\n===== transmit field maps =====\n');
outdir = fullfile(P.dataDir, 'field');
if ~exist(outdir, 'dir'), mkdir(outdir); end

x = P.fieldMapX;  z = P.fieldMapZ;                       % 1 mm 网格
[Xg, Zg] = meshgrid(x, z);
npts = numel(Xg);
pts  = [Xg(:), zeros(npts,1), Zg(:)];        % calc_hp 要求 N x 3 的点坐标

Th = f2_setup(P);
peakAll = zeros(numel(z), numel(x), P.nAngles, 'single');

for i = 1:P.nAngles
    a = P.angles(i);
    d = P.elementX * sin(a) / P.soundSpeed;                 % 与采集时完全相同的发射延时
    xdc_center_focus(Th, [0 0 0]);
    xdc_focus_times(Th, 0, d);
    xdc_apodization(Th, 0, ones(1, P.nElements));    % 时间线 0 s 起生效（不是通道编号）

    % 每个点只关心“峰值有多大”，不关心它出现在哪个时刻：记录长度与时间基准
    % 都不需要对齐，取 max(abs()) 即可。
    hp = calc_hp(Th, pts);                   % [Nt x npts] 声压
    peakAll(:,:,i) = single(reshape(max(abs(hp), [], 1), numel(z), numel(x)));
end

S = struct('x', x, 'z', z, 'peakAll', peakAll);
save(fullfile(outdir, 'maps.mat'), '-struct', 'S', '-v7.3');   % 逐角度声场，用 v7.3
fprintf('saved field maps (%d angles, %d x %d grid)\n', P.nAngles, numel(x), numel(z));
end
