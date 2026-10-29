function T = result_stamp(P, extra)
%RESULT_STAMP  “计算型产物”的输入指纹：声场图、逐点 PSF 图这类确定性结果。
%
%   T = RESULT_STAMP(P) 返回一个结构体，把“重算这份产物会用到”的输入全部装进去。
%   T = RESULT_STAMP(P, EXTRA) 再并进 EXTRA 的字段（网格、细网格步长之类不在 P 里
%   的输入）。EXTRA 不能含 version / stamp 这两个保留名。
%
%   为什么需要它
%   ------------
%   data_stamp 管的是“采集结果”（data/pw、data/foc、data/noise）是否过期；这里的
%   产物不是采集，而是算出来的：stage_fields 用 calc_hp 算 11 个角度的发射声压
%   （约 4 s），run_psfmap 在 165 个点上各跑 11 次发射的完整链路（21–35 s）。两者
%   都只依赖参数，同样的参数必然给出逐位相同的结果，所以重算一遍纯属浪费；但如果
%   参数变了又必须重算，否则图与参数对不上。指纹就是用来分辨这两种情况的。
%
%   取值为什么直接用 data_stamp
%   ---------------------------
%   PSF 图跑的是完整链路（f2_setup → acquire_pw → das_pw → 复合），采集链路那 35
%   项参数一个都不能漏，所以直接复用 data_stamp 的那一套最省事、也不会漏。声场图
%   实际只用到阵列几何、声速、脉冲与角度，取这个超集只会多算、不会拿旧值冒新值。
%   实践中超集并不带来额外的重算：这 35 项里任何一项变了，stamp_guard 会在进入
%   任何阶段之前就停下，根本走不到这里。
%
%   T.version
%   ---------
%   与 data_stamp 的 version 分开：那个管采集算法，这个管“计算产物”的算法。改了
%   stage_fields / psf_map 里影响结果的逻辑（例如加衰减、改窗），把它 +1，旧文件
%   就会被判为过期。
%
%   不跟踪的东西
%   ------------
%   Field II / MATLAB 版本、机器：与 data_stamp 一样只看参数。换了版本想强制重算，
%   用 run_refresh('force')。

T = struct();
T.version = 1;              % 计算方法改动时 +1（见上面 T.version 一节）
T.stamp   = data_stamp(P);  % 采集链路的 35 项参数 + 采集算法版本号

if nargin > 1 && ~isempty(extra)
    f = fieldnames(extra);
    for k = 1:numel(f)
        if any(strcmp(f{k}, {'version', 'stamp'}))
            error('result_stamp:reserved', ...
                'EXTRA 里不能用保留字段名 "%s"（它属于指纹本身）', f{k});
        end
        T.(f{k}) = extra.(f{k});
    end
end
end
