function tf = force_flag(x)
%FORCE_FLAG  把 run_refresh / run_psfmap 的可选参数统一成逻辑量。
%
%   接受三种写法
%     []        或省略 → false（能用缓存就用缓存）
%     'force'        → true （忽略缓存，强制重算）
%     true / 1       → true
%
%   为什么允许字符串：命令行里 run_refresh('force') 比 run_refresh(true) 更容易看懂，
%   也不会被误读成“跑第 1 个阶段”。

if nargin < 1 || isempty(x)
    tf = false;
elseif ischar(x) || (isstring(x) && isscalar(x))
    tf = strcmpi(char(x), 'force');
else
    tf = logical(x);
end
end
