function make_fig_coverage(root)
%MAKE_FIG_COVERAGE  图 F14：11 束平面波照亮的范围，与整幅成像范围的关系。
%
%   MAKE_FIG_COVERAGE(ROOT) 读 data/field/maps.mat（11 个角度在 0.5 mm 网格上的
%   发射峰值声压，每个格点已做 5 个亚样本的平均，见附录 B22），画两幅：
%     (a) 几何：顶部画换能器阵列（与图 F1 同一画法），下面每束平面波的两条边界线、
%         11 束全部照到的区域，以及“已验证均匀”的成像区域（暗红实线框）。这一幅
%         就是整幅模型的示意图：换能器 → 11 束 → 覆盖区 → 成像区域；
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
%   为什么是"幅度相加"而不是"相干相加"
%   --------------------------------
%   DAS 的延时已经把每个角度到该像素的传播相位补齐了，11 个角度的复信号在这一点
%   上同相，所以矢量和的模就等于各模之和（sum |p_alpha|）。换句话说，E 就是图像
%   在那个位置"看到"的发射亮度。
%   反过来，如果相位没有逐像素补齐（例如只按轴上一个固定延时去平移各角度的声场），
%   叠加就会变成干涉条纹——那正是被废弃的旧 F5 所犯的错误，见附录 B13。这一幅
%   与 F5 是同一件事的两种看法：F14 看的是"发射加权平不平"，F5 看的是"成像之后
%   平不平"。
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

% 成像区域框的颜色：场图是 parula（深蓝→青→绿→黄），白框落在黄色背景上几乎看不见，
% 所以两个子图统一改用深红实线——它不在 parula 的色域里，又和 (b) 里黑色的 −1 dB
% 等值线、白色的几何边界都不冲突。
roiCol = [0.72 0.05 0.05];

% 两个子图共用同一套坐标范围，且都用 axis image 锁成 1:1——这样 (a) 里 10° 的边界
% 看上去就是 10°，(b) 的白色边界线与 (a) 的灰区边界也严丝合缝地对齐。
axX = [-24 24];
axZ = [-5.2 60];        % 上方留 5.2 mm 给换能器示意图

fig = figure('Position', [60 60 1250 800], 'Color', 'w');
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

% 换能器模型（与图 F1 同一画法），这样 F14 自己就是整幅模型的示意图：
% 换能器 → 11 束 → 覆盖区 → 成像区域。
draw_array(P, -3.2, 1.1);
text(0, -4.35, sprintf('%d 阵元，间距 %.2f mm，孔径 %.1f mm', ...
     P.nElements, P.elemPitch*1e3, P.arrayLength*1e3), ...
     'HorizontalAlignment', 'center', 'FontSize', 7.5);

plot(roiX, roiZ, '-', 'Color', roiCol, 'LineWidth', 1.8);
text(-9.4, 47.5, '已验证均匀区', 'FontSize', 8);
text(-23, 36, sprintf('11 束全覆盖：\n|x| \\leq %.1f mm\n（z = 50 mm 处）', ...
     half - 50*tan(thMax)), 'FontSize', 8);
xlabel('x (mm)'); ylabel('z (mm)');
axis image;
set(gca, 'XLim', axX, 'YLim', axZ, 'YTick', 0:10:60, 'FontSize', 9, ...
         'Layer', 'top', 'YDir', 'reverse');   % 换能器在上：深度由上至下增大（全项目统一）
title('(a) 11 束平面波的几何覆盖', 'FontSize', 10);

% ---------- (b) 成像加权：11 个角度的幅度和 ----------
nexttile;
E   = sum(double(M.peakAll), 3);             % sum |p_alpha|，DAS 之后图像看到的加权
% 每个 |p_alpha| 都随深度衰减，直接画 E 就只剩一条深度渐变，看不见横向的覆盖。
% 所以按深度归一化到轴线上（E(0,z)）：这个比值就是“该点被几束照到”的度量，
% 11 束全照到时为 1（0 dB），只剩一束时约 1/11（-20.8 dB）。
E0  = E(:, 1+ (numel(xg)-1)/2);              % x = 0 那一列（网格含 0，正中）
EdB = 20*log10(max(E, eps) ./ E0);
[xu, zu, Eu] = zoom3(xg, zg, EdB);
imagesc(xu, zu, Eu); axis image; set(gca, 'YDir', 'reverse'); hold on;  % 换能器在上
colormap(gca, parula); cb = colorbar; cb.Label.String = 'dB（相对同深度轴上）';
caxis([-6 0]);                               % 只看 0..-6 dB：更外面已经离开成像范围
contour(xu, zu, Eu, [-1 -1], 'k--', 'LineWidth', 1.0);   % -1 dB 边界：11 束基本都在
plot(-xAll, zGeo, 'w-', 'LineWidth', 1.2);   % 与 (a) 同一条边界
plot( xAll, zGeo, 'w-', 'LineWidth', 1.2);
plot(roiX, roiZ, '-', 'Color', roiCol, 'LineWidth', 2.0);   % 成像区域：深红实线
text(-17.5, 8.2, '黑虚线：−1 dB 边界；暗红框：成像区域', 'FontSize', 7);
draw_array(P, -3.2, 1.1);                    % 与 (a) 同一个换能器示意，两幅对齐
xlabel('x (mm)'); ylabel('z (mm)');
set(gca, 'XLim', axX, 'YLim', axZ, 'YTick', 0:10:60, 'FontSize', 9);
title('(b) 成像加权 \Sigma_\alpha |p_\alpha|（0 dB = 11 束全覆盖）', 'FontSize', 10);

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
function draw_array(P, zBottom, h)
%DRAW_ARRAY  换能器示意图：每隔 4 个阵元画一个小矩形。
%
%   128 个全画出来就是一条黑杠（0.33 mm 间距在 40 mm 宽的画面上不足 1 像素），
%   所以与图 F1 用同一画法。zBottom 是矩形下沿的深度（mm），h 是高度（mm）；
%   调用前坐标轴必须处于 hold on 状态。

for k = 1:4:P.nElements
    rectangle('Position', [(P.elementX(k)-P.elemWidth/2)*1e3, zBottom, ...
                           P.elemWidth*1e3, h], ...
              'FaceColor', [0.35 0.35 0.35], 'EdgeColor', 'none');
end
end

% =========================================================================
function [xv, zv, Vu] = zoom3(xg, zg, V)
%ZOOM3  3 倍双线性插值，只为显示（0.5 mm 网格直接放大成 0.17 mm 的显示像素）。
xv = linspace(xg(1), xg(end), 3*numel(xg));
zv = linspace(zg(1), zg(end), 3*numel(zg));
[XU, ZU] = meshgrid(xv, zv);
Vu = interp2(xg, zg, V, XU, ZU, 'linear');
end
