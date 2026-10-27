function make_fig_coverage(root)
%MAKE_FIG_COVERAGE  图 F14：11 束平面波照亮的范围，与整幅成像范围的关系。
%
%   MAKE_FIG_COVERAGE(ROOT) 读 data/field/maps.mat（11 个角度在 1 mm 网格上的
%   发射峰值声压），画两幅：
%     (a) 几何：每束平面波的两条边界线、11 束全部照到的区域、孔径、以及
%         “已验证均匀”的成像区域（虚线框）；
%     (b) 成像加权 E = sum_alpha |p_alpha|：DAS 把每个角度的传播相位逐像素补掉
%         之后，图像实际看到的发射加权（幅度相加，不是相干相加）。
%
%   为什么要这一幅
%   --------------
%   11 束的公共照明区随深度收窄：孔径半宽 21.12 mm、最大偏转 10°，在深度 z 处
%   公共区是 |x| <= 21.12 - z*tan(10°) mm，z=60 mm 时只剩 ±10.5 mm。成像区域
%   （|x| <= 10 mm、z = 15..50 mm）就落在这块公共区里，所以那里的加权是平的
%   ——这正是 F5 逐点验证幅值均匀的前提。
%
%   用法：make_fig_coverage

if nargin < 1 || isempty(root)
    root = fileparts(fileparts(mfilename('fullpath')));
end
P = params(root);
M = load(fullfile(P.dataDir, 'field', 'maps.mat'));
F = fullfile(root, 'figures');
if ~exist(F, 'dir'), mkdir(F); end

xg = M.x*1e3;  zg = M.z*1e3;                 % mm
half  = P.arrayLength/2 * 1e3;                         % 孔径半宽 21.12 mm
thMax = max(abs(P.angles));                  % 10°
roiX  = [-10 10 10 -10 -10];                 % 已验证均匀的成像区域（与 F5 同一块）
roiZ  = [15 15 50 50 15];

% 11 束的公共照明区半宽（几何）：x_all(z) = half - z*tan(10°)
zGeo  = linspace(0, 60, 200);
xAll  = half - zGeo*tan(thMax);

fig = figure('Position', [60 60 1300 700], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---------- (a) 几何：谁照到哪儿 ----------
nexttile; hold on; box on;
cmap = parula(numel(P.angles));
for q = 1:numel(P.angles)
    th = P.angles(q);
    plot(-half + zGeo*tan(th), zGeo, '-', 'Color', cmap(q,:), 'LineWidth', 0.9);
    plot( half + zGeo*tan(th), zGeo, '-', 'Color', cmap(q,:), 'LineWidth', 0.9);
end
patch([-xAll, fliplr(xAll)], [zGeo, fliplr(zGeo)], [0.25 0.25 0.25], ...
      'FaceAlpha', 0.12, 'EdgeColor', 'none');            % 11 束全部照到
patch([-half half half -half], [-0.9 -0.9 0 0], [0.6 0.6 0.6], 'EdgeColor', 'none');
plot(roiX, roiZ, 'k--', 'LineWidth', 1.5);
text(-9.4, 47.5, 'verified region', 'FontSize', 8);
text(-23, 36, sprintf('all 11 beams:\n|x| \\leq %.1f mm\nat z = 50 mm', ...
     half - 50*tan(thMax)), 'FontSize', 8);
text(-23, 12, 'aperture 42.2 mm', 'FontSize', 8);
xlabel('x (mm)'); ylabel('z (mm)');
set(gca, 'XLim', [-24 24], 'YLim', [0 60], 'FontSize', 9, 'Layer', 'top');
title('(a) geometric coverage of the 11 beams', 'FontSize', 10);

% ---------- (b) 成像加权：11 个角度的幅度和 ----------
nexttile;
E   = sum(double(M.peakAll), 3);             % sum |p_alpha|，DAS 之后图像看到的加权
% 每个 |p_alpha| 都随深度衰减，直接画 E 就只剩一条深度渐变，看不见横向的覆盖。
% 所以按深度归一化到轴线上（E(0,z)）：这个比值就是“该点被几束照到”的度量，
% 11 束全照到时为 1（0 dB），只剩一束时约 1/11（-20.8 dB）。
E0  = E(:, 1+ (numel(xg)-1)/2);              % x = 0 那一列（网格含 0，正中）
EdB = 20*log10(max(E, eps) ./ E0);
[xu, zu, Eu] = zoom3(xg, zg, EdB);
imagesc(xu, zu, Eu); axis xy image; hold on;
colormap(gca, parula); cb = colorbar; cb.Label.String = 'dB re on-axis, same depth';
caxis([-6 0]);                               % 只看 0..-6 dB：更外面已经离开成像范围
contour(xu, zu, Eu, [-1 -1], 'k--', 'LineWidth', 1.0);   % -1 dB 边界：11 束基本都在
plot(-xAll, zGeo, 'w-', 'LineWidth', 1.2);   % 与 (a) 同一条边界
plot( xAll, zGeo, 'w-', 'LineWidth', 1.2);
plot(roiX, roiZ, 'w--', 'LineWidth', 1.5);
text(-19.5, 6.5, 'dashed: -1 dB boundary', 'FontSize', 7.5);
xlabel('x (mm)'); ylabel('z (mm)'); set(gca, 'FontSize', 9);
title('(b) imaging weighting \Sigma_\alpha |p_\alpha| (0 dB = all 11 beams)', 'FontSize', 10);

% 框内起伏：把 ROI 里归一化后的最大/最小差打印出来，报告直接引用
[XM, ZM] = meshgrid(xg, zg);
in = abs(XM) <= 10 & ZM >= 15 & ZM <= 50;
fprintf('ROI (|x|<=10, z=15..50): E/E0 spread = %.2f dB (min %.2f dB)\n', ...
        max(EdB(in)) - min(EdB(in)), min(EdB(in)));

out = fullfile(F, 'F14_field_coverage.png');
save_fig(fig, out); close(fig);
fprintf('F14 saved: %s\n', out);
end

% =========================================================================
function [xv, zv, Vu] = zoom3(xg, zg, V)
%ZOOM3  3 倍双线性插值，只为显示（1 mm 网格直接放大成块状纹理）。
xv = linspace(xg(1), xg(end), 3*numel(xg));
zv = linspace(zg(1), zg(end), 3*numel(zg));
[XU, ZU] = meshgrid(xv, zv);
Vu = interp2(xg, zg, V, XU, ZU, 'linear');
end
