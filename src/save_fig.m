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

if nargin < 3 || isempty(dpi)
    dpi = 200;                              % 默认值：嵌入 Word 时清晰度与体积的折中
end

d = fileparts(path);                        % 确保目标目录存在
if ~isempty(d) && ~exist(d, 'dir')
    mkdir(d);
end

try
    exportgraphics(fig, path, 'Resolution', dpi);   % 首选：保持屏幕版面
catch
    print(fig, path, '-dpng', sprintf('-r%d', dpi));% 兜底：老式导出路径
end
end
