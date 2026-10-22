function save_fig(fig, path, dpi)
%SAVE_FIG  Write a figure as PNG (and keep the handle open for the caller).

if nargin < 3 || isempty(dpi), dpi = 200; end
d = fileparts(path);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
try
    exportgraphics(fig, path, 'Resolution', dpi);
catch
    print(fig, path, '-dpng', sprintf('-r%d', dpi));
end
end
