function lines = stamp_diff(a, b, prefix)
%STAMP_DIFF  逐字段比较两个参数结构，返回人类可读的差异列表。
%
%   LINES = STAMP_DIFF(A, B) 比较结构体 A（旧）与 B（新），返回字符串元胞数组，
%   每个元素形如 "pitch: 0.00033 -> 0.00034"。没有差异时返回空元胞。
%
%   LINES = STAMP_DIFF(A, B, PREFIX) 给每条差异加前缀，递归时用。

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
if isnumeric(a) && isnumeric(b) && isequal(size(a), size(b)) && numel(a) > 8
    d = max(abs(double(a(:)) - double(b(:))));
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
