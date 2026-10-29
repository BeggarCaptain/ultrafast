function save_fig(fig, path, dpi)
%SAVE_FIG  把图窗按指定分辨率保存为 PNG 文件。
%
%   SAVE_FIG(FIG, PATH, DPI) 将图窗句柄 FIG 导出到 PATH。本函数不关闭图窗：
%   句柄仍归调用者所有，由调用者决定何时 close（有的调用者保存后还要继续加
%   子图，提前关闭会出错）。
%
%   输入
%     fig  ：图窗句柄
%     path ：目标文件名（父目录不存在会自动创建）
%     dpi  ：栅格分辨率，默认 200 dpi（要印刷质量可用 300）
%
%   优先用 exportgraphics：它按屏幕上的实际版面导出，不会重新排版、不会把
%   坐标轴标签裁掉。print 只是兜底——极少数图（例如超大 uipanel）在新版本上
%   仍可能只支持 print。
%
%   导出前统一换成含汉字字形的字体
%   ----------------------------------
%   MATLAB 默认的 Helvetica 没有汉字，中文标注会画成一排方框。项目里所有图都
%   经本函数导出，所以字体只在这里设一次；各张图只管写中文，不必重复指定字体。
%   注意字体要设到图内**所有**文字对象（标题、刻度、图例、色标标签），只设 gca
%   的 FontName 管不到色标与图例。
%
%   为什么默认 200 dpi
%   -----------------
%   报告里的图宽约 16 cm，200 dpi 下约 1260 px，屏幕阅读与打印都够清晰，单图
%   通常不到 1 MB；改 300 dpi 会让 14 张图的总体积与生成时间都接近翻倍，而这个
%   项目不投稿印刷版。要印刷时在调用处传 300 即可。
%
%   为什么要补白边
%   -------------
%   exportgraphics 会把画布裁到内容边界，四周只剩 2 px：标题与刻度数字几乎贴着
%   边，插进 Word 后容易被表格线或页边吃掉。R2024a 的 exportgraphics 没有
%   Padding 选项（传了会报“非法的选项”），imwrite 也没有补边参数，所以这里直接
%   复制边缘行列手工补——只增加白边，不改动任何一个原始像素。

if nargin < 3 || isempty(dpi)
    dpi = 200;                              % 默认值：嵌入 Word 时清晰度与体积的折中
end

% 中文字体：按优先级取第一个系统里装了的（Windows 自带雅黑，缺失时退回思源/黑体）
prefer = {'Microsoft YaHei UI', 'Microsoft YaHei', 'Noto Sans SC', 'SimHei', 'SimSun'};
have   = listfonts;
fontCN = prefer{1};
for k = 1:numel(prefer)
    if any(strcmpi(have, prefer{k}))
        fontCN = prefer{k};
        break;
    end
end
set(findall(fig, '-property', 'FontName'), 'FontName', fontCN);   % 含标题、图例、色标

d = fileparts(path);                        % 确保目标目录存在
if ~isempty(d) && ~exist(d, 'dir')
    mkdir(d);
end

try
    exportgraphics(fig, path, 'Resolution', dpi);   % 首选：保持屏幕版面
catch
    print(fig, path, '-dpng', sprintf('-r%d', dpi));% 兜底：老式导出路径
end

% exportgraphics 会把画布裁到内容边界，四周只剩 2 px：标题与刻度数字几乎贴着边，
% 插进 Word 后容易被表格线或页边吃掉。这里统一补一圈白边（按 dpi 缩放）。
I = imread(path);
if ismatrix(I), I = repmat(I, 1, 1, 3); end        % 灰度 PNG 也补成三通道
pad = max(6, round(dpi/20));
I = I([ones(1,pad) 1:size(I,1) size(I,1)*ones(1,pad)], ...
      [ones(1,pad) 1:size(I,2) size(I,2)*ones(1,pad)], :);
imwrite(I, path);
end
