function make_fig_pointtargets()
%MAKE_FIG_POINTTARGETS  图 F14：8 个点目标的三方法成像 + 真实位置标注。
%
%   直接读取已经采集好的 P3 数据（data/pw/P3_a###.mat、data/foc/P3_l###.mat），
%   不重新跑任何 Field II 仿真。做三件事：
%     1. 把 PW（0°）、C11（11 个角度复合）、单焦点 三种图像画出来；
%     2. 在每个目标的真实坐标上叠加标记（红圈/十字），用眼睛就能判断成像是否
%        落在正确位置；
%     3. 逐目标搜索包络峰值，输出“实测 − 真实”的定位误差与横向/轴向宽度，
%        同时写 data/res/table_pointtargets.csv。
%
%   为什么用点目标而不是均匀散斑来看“成像准不准”
%   ------------------------------------------------
%   均匀散斑体模的图像本来就是一格一格的斑点（相干干涉的结果），看起来像噪声
%   是正确物理，不是仿真出错；它能回答的是分辨率/对比度，不能回答“位置准不准”。
%   点目标的真实坐标是已知的，所以适合用来核对几何定位。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);

% ---------------------------------------------------------------- 读数据
% PW：只取 0° 那一次
S0  = load(fullfile(P.dataDir, 'pw', sprintf('P3_a%03d.mat', P.iZero)));
rfPW = double(S0.img);

% C11：11 个角度的射频直接相加（仍在射频域）
rfC = zeros(size(rfPW), 'single');
for i = 1:P.nAngles
    S = load(fullfile(P.dataDir, 'pw', sprintf('P3_a%03d.mat', i)));
    rfC = rfC + S.img;
end
rfC = double(rfC);

% 单焦点：256 条线拼成图像
rfF = zeros(numel(P.z), numel(P.xpsf), 'single');
for k = 1:numel(P.xpsf)
    S = load(fullfile(P.dataDir, 'foc', sprintf('P3_l%03d.mat', k)));
    rfF(:,k) = S.img;
end
rfF = double(rfF);

envPW = abs(hilbert(rfPW));
envC  = abs(hilbert(rfC));
envF  = abs(hilbert(rfF));

% ------------------------------------------------------------ 真实目标
[pos, ~] = phantom_points(P);
xt = pos(:,1);  zt = pos(:,3);
nt = numel(xt);

% ------------------------------------------- 逐目标的定位误差与宽度
% 在真实位置附近开一个窗口，找包络最大值作为“实测位置”。窗口 3x2 mm 足够
% 容纳主瓣，又不会框进相邻目标（最近的两个目标相距 10 mm）。
%
% 两种读数都记录：
%   raw —— 直接取峰值所在像素，误差被网格量化到 dx = 0.171 mm 的整数倍
%          （x = 0 正好落在两列之间，所以 raw 读数会稳定偏 ±dx/2 = ±86 µm）；
%   fit —— 在峰值两侧各取一点做三点抛物线插值，得到亚像素位置。
% 报告里用 fit，raw 用来说明“误差来自网格，不是波束形成”。
res = zeros(nt, 7);                       % [dxRaw dzRaw dxFit dzFit w10 w6]
for q = 1:nt
    ix = find(abs(P.x - xt(q)) <= 1.5e-3);
    iz = find(abs(P.z - zt(q)) <= 1.0e-3);
    if isempty(ix) || isempty(iz), continue; end
    E = envC(iz, ix);
    [~, im] = max(E(:));
    [r, c] = ind2sub(size(E), im);

    % ---- 粗读数：像素中心
    xRaw = P.x(ix(c));  zRaw = P.z(iz(r));

    % ---- 亚像素：三点抛物线（只有峰值不在窗口边缘时才能做）
    r0 = iz(1) + r - 1;  c0 = ix(1) + c - 1;
    xFit = P.x(c0);  zFit = P.z(r0);
    if c0 > 1 && c0 < numel(P.x)
        a = envC(r0, c0-1);  b = envC(r0, c0);  d = envC(r0, c0+1);
        den = a - 2*b + d;
        if abs(den) > 0, xFit = P.x(c0) + 0.5*(a-d)/den*P.dx; end
    end
    if r0 > 1 && r0 < numel(P.z)
        a = envC(r0-1, c0);  b = envC(r0, c0);  d = envC(r0+1, c0);
        den = a - 2*b + d;
        if abs(den) > 0, zFit = P.z(r0) + 0.5*(a-d)/den*P.dz; end
    end

    res(q,1) = xRaw - xt(q);  res(q,2) = zRaw - zt(q);
    res(q,3) = xFit - xt(q);  res(q,4) = zFit - zt(q);
end

% 轴向/横向宽度用统一的指标函数（与报告 T-B/T-C 同口径）
for q = 1:nt
    m = metrics_resolution(single(envC), P.x, P.z, xt(q), zt(q));
    res(q,5) = m.w10;  res(q,6) = m.w6;
end

fprintf('\n point-target localisation (C11)\n');
fprintf('   #    true (mm)      raw peak (mm)   err_raw  err_fit  w10    w6\n');
for q = 1:nt
    fprintf('  %2d  (%6.1f,%5.1f)  (%6.2f,%6.2f)  (%+4.0f,%+3.0f)  (%+5.0f,%+4.0f)  %.2f  %.2f\n', ...
        q, xt(q)*1e3, zt(q)*1e3, ...
        xt(q)*1e3 + res(q,1)*1e3, zt(q)*1e3 + res(q,2)*1e3, ...
        res(q,1)*1e6, res(q,2)*1e6, res(q,3)*1e6, res(q,4)*1e6, ...
        res(q,5)*1e3, res(q,6)*1e3);
end
fprintf('  raw : |dx| max %4.0f um, |dz| max %3.0f um  (grid dx/2 = %.0f um)\n', ...
    max(abs(res(:,1)))*1e6, max(abs(res(:,2)))*1e6, P.dx/2*1e6);
fprintf('  fit : |dx| max %4.0f um, |dz| max %3.0f um\n', ...
    max(abs(res(:,3)))*1e6, max(abs(res(:,4)))*1e6);

% CSV
fid = fopen(fullfile(P.dataDir,'res','table_pointtargets.csv'), 'w');
fprintf(fid, 'index,x_true_mm,z_true_mm,x_fit_mm,z_fit_mm,dx_raw_um,dz_raw_um,dx_fit_um,dz_fit_um,w10_mm,w6_mm\n');
for q = 1:nt
    fprintf(fid, '%d,%.2f,%.2f,%.2f,%.2f,%.0f,%.0f,%.0f,%.0f,%.3f,%.3f\n', q, ...
        xt(q)*1e3, zt(q)*1e3, xt(q)*1e3+res(q,3)*1e3, zt(q)*1e3+res(q,4)*1e3, ...
        res(q,1)*1e6, res(q,2)*1e6, res(q,3)*1e6, res(q,4)*1e6, ...
        res(q,5)*1e3, res(q,6)*1e3);
end
fclose(fid);

% ------------------------------------------------------------------ 画图
dbPW = to_db(envPW, P.drDB);
dbC  = to_db(envC,  P.drDB);
dbF  = to_db(envF,  P.drDB);

fig = figure('Position', [80 80 1150 760]);

% ---- 上排：全场 ----
imgs = {dbPW, dbC, dbF};
ttls = {'PW (n=1)', 'C11 (n=11)', 'focused (256 lines)'};
xs   = {P.x, P.x, P.xpsf};
for q = 1:3
    subplot(2,3,q);
    imagesc(xs{q}*1e3, P.z*1e3, imgs{q}); axis image; colormap(gca, gray);
    caxis([-P.drDB 0]); hold on;
    plot(xt*1e3, zt*1e3, 'ro', 'MarkerSize', 12, 'LineWidth', 1.2);   % 真实位置
    xlabel('x (mm)'); ylabel('z (mm)'); title(ttls{q});
end

% ---- 下排：30 mm 处放大，十字标真实位置 ----
z0 = 30e-3;  x0 = 0;
iz = abs(P.z - z0) <= 2.5e-3;
ix = abs(P.x - x0) <= 3e-3;
for q = 1:3
    subplot(2,3,3+q);
    imagesc(xs{q}(ix)*1e3, P.z(iz)*1e3, imgs{q}(iz, ix)); axis image;
    colormap(gca, gray); caxis([-40 0]); hold on;
    plot(x0*1e3, z0*1e3, 'r+', 'MarkerSize', 14, 'LineWidth', 1.4);
    xlabel('x (mm)'); ylabel('z (mm)');
    title(sprintf('%s, zoom at (0, 30) mm', ttls{q}), 'FontSize', 9);
end
save_fig(fig, fullfile(P.figDir, 'F14_point_targets.png')); close(fig);
fprintf('saved %s\n', fullfile(P.figDir, 'F14_point_targets.png'));
end
