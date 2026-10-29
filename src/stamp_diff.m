function lines = stamp_diff(a, b, prefix)
%STAMP_DIFF  逐字段比较两个参数结构，返回人类可读的差异列表。
%
%   LINES = STAMP_DIFF(A, B) 比较结构体 A（旧）与 B（新），返回字符串元胞数组，
%   每个元素形如 "pitch: 0.00033 -> 0.00034"。没有差异时返回空元胞。
%
%   LINES = STAMP_DIFF(A, B, PREFIX) 给每条差异加前缀，递归时用。
%
%   为什么不用 isequal 一句话解决
%   ----------------------------
%   isequal 只能回答"一样不一样"，而人需要的是"哪个参数变了、从多少变到多少"。
%   stamp_guard 撞到不一致时要立刻让人判断"这是误改还是我确实想重跑"，所以差异
%   必须具体到字段名与数值。
%
%   怎么比较
%   --------
%   递归展开结构体（指纹是 S.stamp.fields.<参数名> 这种嵌套结构），用点号拼出
%   完整路径，例如 "stamp.fields.angles"。每个字段三种情况：
%     * 旧有新无 → "新增（新值 …）"；
%     * 新有旧无 → "已删除（旧值 …）"；
%     * 都有     → 值不同就报出来。
%
%   为什么长数组只报尺寸与 max|差值|
%   ---------------------------------
%   角度表是 11 个数、点目标坐标是 8 个——都在一屏之内，可以直接列。像素网格、
%   成像坐标这类几百上千个数的向量列出来没法看，而"同尺寸 [1 11]，max|差值|
%   = 0.175"已经足够判断是否实质变化（0.175 显然是有人动了角度，1e-16 则只是
%   浮点噪声）。
%
%   为什么用 isequaln 而不是 isequal
%   --------------------------------
%   NaN 在本项目里是合法取值（例如 0° 的倾角比值无定义，存的就是 NaN）。
%   isequal(NaN, NaN) 为假，会把"两个都是 NaN"误报成"变了"；isequaln 把 NaN
%   当作相等。

if nargin < 3 || isempty(prefix)
    prefix = '';
end
lines = {};

% 非结构体（标量、数组、字符串）：直接比较
if ~isstruct(a) || ~isstruct(b)
    if ~isequaln(a, b)
        lines{end+1} = sprintf('%s: %s -> %s', prefix, stamp_short(a), stamp_short(b));
    end
    return;
end

fa = fieldnames(a);
fb = fieldnames(b);
for k = 1:numel(fa)
    n = fa{k};
    if isempty(prefix)
        p = n;
    else
        p = [prefix '.' n];
    end
    if ~isfield(b, n)
        lines{end+1} = sprintf('%s: 已删除（旧值 %s）', p, stamp_short(a.(n))); %#ok<AGROW>
    elseif isstruct(a.(n)) && isstruct(b.(n))
        lines = [lines, stamp_diff(a.(n), b.(n), p)]; %#ok<AGROW>
    elseif ~isequaln(a.(n), b.(n))
        lines{end+1} = sprintf('%s: %s', p, pair_str(a.(n), b.(n))); %#ok<AGROW>
    end
end
for k = 1:numel(fb)
    if ~isfield(a, fb{k})
        lines{end+1} = sprintf('%s.%s: 新增（新值 %s）', prefix, fb{k}, ...
            stamp_short(b.(fb{k}))); %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function s = pair_str(a, b)
%PAIR_STR  一对值的差异描述：长数组只报尺寸与最大绝对差，短值直接列出来。
%   "长"的门槛取 8 个元素：参数里最长的短向量是 11 个角度、8 个点目标坐标，
%   它们恰好能一行列完；再长就换成长度摘要，否则提示信息会被数字淹没。
if isnumeric(a) && isnumeric(b) && isequal(size(a), size(b)) && numel(a) > 8
    d = max(abs(double(a(:)) - double(b(:))));
    % 尺寸不同（例如角度从 11 个改成 9 个）会落到 else 分支，直接列出完整向量
    % ——那种情况下"哪个角度被删了"比"最大差值多少"更重要。
    s = sprintf('同尺寸 %s，max|差值| = %.3g', mat2str(size(a)), d);
else
    s = sprintf('%s -> %s', stamp_short(a), stamp_short(b));
end
end

% -------------------------------------------------------------------------
function s = stamp_short(x)
%STAMP_SHORT  把任意值压成一行短字符串，用于差异提示。
if ischar(x)
    s = x;
elseif isstring(x)
    s = char(x);
elseif isscalar(x) && isnumeric(x)
    s = sprintf('%.6g', x);
elseif isnumeric(x) && numel(x) <= 8
    s = mat2str(x, 6);
elseif isnumeric(x)
    s = sprintf('[%s 数值]', strjoin(cellfun(@num2str, num2cell(size(x)), ...
        'UniformOutput', false), 'x'));
else
    s = sprintf('<%s>', class(x));
end
end
