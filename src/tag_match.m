function ok = tag_match(file, tag)
%TAG_MATCH  磁盘上的产物是不是“当前参数算出来的那一份”。
%
%   OK = TAG_MATCH(FILE, TAG) 读 FILE 里存的指纹并与 TAG 逐项比较。文件不存在、
%   文件里没有 tag 字段（旧版本写的）、或任一项不同，都返回 false——宁可重算一次，
%   也不用可能是旧参数的结果。
%
%   TAG 由 src/result_stamp.m 生成。逐项比较用 isequal，所以参数是精确相等才算
%   命中：这些参数都由 config/params.m 里的字面量和整数索引确定性地算出来，同一
%   套配置每次得到的双精度值完全相同，不存在“差不多相等”的情况。
%
%   只为读 tag 就把整个文件读进来：stage_fields 的 maps.mat 约 0.1 MB，psfmap.mat
%   更小，代价可以忽略；这样也不必为 v7.3 文件单独做部分读取。

ok = false;
if exist(file, 'file') ~= 2
    return;
end
S = load(file);
ok = isfield(S, 'tag') && isequal(S.tag, tag);
end
