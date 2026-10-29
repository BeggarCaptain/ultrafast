function make_fig_pointtargets()
%MAKE_FIG_POINTTARGETS  图 F13：8 个点目标的三方法成像 + 真实位置标注。
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
%
%   这张图验证的是什么
%   ----------------
%   它验证的不是分辨率，而是**几何**：延时模型（发射 + 接收 + Field II 时间轴
%   修正 + dt）与 DAS 加起来，是否把目标放回了它真实所在的位置。如果某项符号写
%   反、或漏掉 +x1 sin(alpha)/c，图像会整体偏移或随深度剪切——在这种“已知真值”
%   的图上立刻就能看出来，而在散斑图上看不出来。
%
%   为什么误差要报两个数（raw 与 fit）
%   --------------------------------
%   raw 是“峰值落在哪个像素”，它的误差被网格量化：横向像素 0.171 mm，x = 0
%   正好落在两列之间，所以 raw 会稳定偏 ±86 µm。那是网格的刻度，不是算法的精度。
%   fit 在峰值两侧各取一点做三点抛物线插值（主瓣顶部二阶展开就是抛物线），能
%   给出亚像素位置，实测误差 3–12 µm。报告用 fit，raw 留在表里用来说明“误差
%   主要来自网格”。
%
%   为什么 C11 与常规聚焦各报一组
%   ---------------------------
%   C11 用 P.imageX、聚焦用 P.psfLineX（同为 256 线、0.171 mm，但轴向量不同，
%   不能混用）。图 F13 的聚焦面板正是暴露 das_focused 缺陷的地方：发射路程原写
%   成 z、少了 |F|-zf，使离轴线的整幅图偏深 (|F|-zf)/2（|x|=12 mm 处 1.16 mm，
%   轴上为 0——只看轴上六个点永远发现不了）。修正后两组误差都在 0.09 mm 以内，
%   这组数就是那次修正的验收证据（附录 B20）。
%
%   这里的 C11 与 stage_analysis 里的 C11 是同一算法
%   ------------------------------------------------
%   都是“11 个角度的射频直接相加，再取包络”，所以 F13 展示的就是最终图像，不是
%   另一套近似。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);

% ---------------------------------------------------------------- 读数据
% PW：只取 0° 那一次
S0  = load(fullfile(P.dataDir, 'pw', sprintf('P3_a%03d.mat', P.zeroAngleIndex)));
rfPW = double(S0.img);

% C11：11 个角度的射频直接相加（仍在射频域）
rfC = zeros(size(rfPW), 'single');
for i = 1:P.nAngles
    S = load(fullfile(P.dataDir, 'pw', sprintf('P3_a%03d.mat', i)));
    rfC = rfC + S.img;
end
rfC = double(rfC);

% 单焦点：256 条线拼成图像
rfF = zeros(numel(P.imageZ), numel(P.psfLineX), 'single');
for k = 1:numel(P.psfLineX)
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
res = loc_err(envC, P.imageX,   P, xt, zt);      % C11（256 线，0.171 mm）
resF = loc_err(envF, P.psfLineX, P, xt, zt);     % 常规聚焦（256 线，0.171 mm）

% 轴向/横向宽度用统一的指标函数（与报告 T-B 同口径）
for q = 1:nt
    m = metrics_resolution(single(envC), P.imageX, P.imageZ, xt(q), zt(q));
    res(q,5) = m.wLat;  res(q,6) = m.wAx;
end

fprintf('\n point-target localisation\n');
fprintf('   #    true (mm)       C11: err_raw   err_fit     聚焦: err_raw   err_fit    wLat  wAx\n');
for q = 1:nt
    fprintf('  %2d  (%6.1f,%5.1f)  (%+4.0f,%+3.0f) (%+5.0f,%+4.0f)   (%+4.0f,%+3.0f) (%+5.0f,%+4.0f)  %.2f  %.2f\n', ...
        q, xt(q)*1e3, zt(q)*1e3, ...
        res(q,1)*1e6,  res(q,2)*1e6,  res(q,3)*1e6,  res(q,4)*1e6, ...
        resF(q,1)*1e6, resF(q,2)*1e6, resF(q,3)*1e6, resF(q,4)*1e6, ...
        res(q,5)*1e3,  res(q,6)*1e3);
end
fprintf('  C11   raw |dx| %.0f um, |dz| %.0f um (grid dx/2 = %.0f um) | fit |dx| %.0f, |dz| %.0f um\n', ...
    max(abs(res(:,1)))*1e6,  max(abs(res(:,2)))*1e6,  P.pixelPitchX/2*1e6, ...
    max(abs(res(:,3)))*1e6,  max(abs(res(:,4)))*1e6);
fprintf('  聚焦  raw |dx| %.0f um, |dz| %.0f um                              | fit |dx| %.0f, |dz| %.0f um\n', ...
    max(abs(resF(:,1)))*1e6, max(abs(resF(:,2)))*1e6, ...
    max(abs(resF(:,3)))*1e6, max(abs(resF(:,4)))*1e6);

% CSV（聚焦的定位误差追加在后面四列，旧的读取方式不受影响）
fid = fopen(fullfile(P.dataDir,'res','table_pointtargets.csv'), 'w');
fprintf(fid, 'index,x_true_mm,z_true_mm,x_fit_mm,z_fit_mm,dx_raw_um,dz_raw_um,dx_fit_um,dz_fit_um,wLat_mm,wAx_mm,x_foc_mm,z_foc_mm,dx_foc_um,dz_foc_um\n');
for q = 1:nt
    fprintf(fid, '%d,%.2f,%.2f,%.2f,%.2f,%.0f,%.0f,%.0f,%.0f,%.3f,%.3f,%.2f,%.2f,%.0f,%.0f\n', q, ...
        xt(q)*1e3, zt(q)*1e3, xt(q)*1e3+res(q,3)*1e3, zt(q)*1e3+res(q,4)*1e3, ...
        res(q,1)*1e6, res(q,2)*1e6, res(q,3)*1e6, res(q,4)*1e6, ...
        res(q,5)*1e3, res(q,6)*1e3, ...
        xt(q)*1e3+resF(q,3)*1e3, zt(q)*1e3+resF(q,4)*1e3, resF(q,3)*1e6, resF(q,4)*1e6);
end
fclose(fid);

% ------------------------------------------------------------------ 画图
dbPW = to_db(envPW, P.displayRangeDB);
dbC  = to_db(envC,  P.displayRangeDB);
dbF  = to_db(envF,  P.displayRangeDB);

fig = figure('Position', [80 80 1150 760]);

% ---- 上排：全场 ----
imgs = {dbPW, dbC, dbF};
ttls = {'单平面波', '11 角度复合', '常规聚焦（256 线）'};
xs   = {P.imageX, P.imageX, P.psfLineX};
for q = 1:3
    subplot(2,3,q);
    imagesc(xs{q}*1e3, P.imageZ*1e3, imgs{q}); axis image; colormap(gca, gray);
    caxis([-P.displayRangeDB 0]); hold on;
    plot(xt*1e3, zt*1e3, 'ro', 'MarkerSize', 12, 'LineWidth', 1.2);   % 真实位置
    xlabel('x (mm)'); ylabel('z (mm)'); title(ttls{q});
end

% ---- 下排：30 mm 处放大，十字标真实位置 ----
z0 = 30e-3;  x0 = 0;
iz = abs(P.imageZ - z0) <= 2.5e-3;
ix = abs(P.imageX - x0) <= 3e-3;
for q = 1:3
    subplot(2,3,3+q);
    imagesc(xs{q}(ix)*1e3, P.imageZ(iz)*1e3, imgs{q}(iz, ix)); axis image;
    colormap(gca, gray); caxis([-40 0]); hold on;
    plot(x0*1e3, z0*1e3, 'r+', 'MarkerSize', 14, 'LineWidth', 1.4);
    xlabel('x (mm)'); ylabel('z (mm)');
    title(sprintf('%s\n（30 mm 处放大）', ttls{q}), 'FontSize', 9);
end
save_fig(fig, fullfile(P.figureDir, 'F13_point_targets.png')); close(fig);
fprintf('saved %s\n', fullfile(P.figureDir, 'F13_point_targets.png'));
end

% =========================================================================
function r = loc_err(env, xv, P, xt, zt)
%LOC_ERR  逐目标的定位误差：[dxRaw dzRaw dxFit dzFit]，单位 m。
%   窗口、读数、亚像素插值与主流程完全一致，只是把图像和它的横轴换掉：
%   C11 用 P.imageX、常规聚焦用 P.psfLineX（两者都是 256 线、0.171 mm，但轴向量
%   不同，不能混用）。之所以两种方法都算，是因为图 F13 的聚焦面板曾暴露出
%   das_focused 的发射项少了一项、整条线偏深 (|F|-zf)/2（附录 B20）——定位误差
%   必须两种方法各报一组，才能证明那一项确实修好了。
r = zeros(numel(xt), 4);
for q = 1:numel(xt)
    ix = find(abs(xv - xt(q)) <= 1.5e-3);
    iz = find(abs(P.imageZ - zt(q)) <= 1.0e-3);
    if isempty(ix) || isempty(iz), continue; end
    E = env(iz, ix);
    [~, im] = max(E(:));
    [rr, cc] = ind2sub(size(E), im);

    xRaw = xv(ix(cc));  zRaw = P.imageZ(iz(rr));          % 粗读数：像素中心
    r0 = iz(1) + rr - 1;  c0 = ix(1) + cc - 1;            % 全图下标
    xFit = xv(c0);  zFit = P.imageZ(r0);
    if c0 > 1 && c0 < numel(xv)
        a = env(r0, c0-1);  b = env(r0, c0);  d = env(r0, c0+1);
        den = a - 2*b + d;
        if abs(den) > 0, xFit = xv(c0) + 0.5*(a-d)/den*(xv(2)-xv(1)); end
    end
    if r0 > 1 && r0 < numel(P.imageZ)
        a = env(r0-1, c0);  b = env(r0, c0);  d = env(r0+1, c0);
        den = a - 2*b + d;
        if abs(den) > 0, zFit = P.imageZ(r0) + 0.5*(a-d)/den*P.pixelPitchZ; end
    end
    r(q,:) = [xRaw - xt(q), zRaw - zt(q), xFit - xt(q), zFit - zt(q)];
end
end
