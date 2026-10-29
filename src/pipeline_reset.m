function files = pipeline_reset(P, dryRun)
%PIPELINE_RESET  冷启动：删除全部由 Field II 采集得到的缓存。
%
%   FILES = PIPELINE_RESET(P)        删除并返回被删文件列表
%   FILES = PIPELINE_RESET(P, true)  只列清单，不动磁盘（dry run）
%
%   删除范围（都是可以重算的中间结果）
%   ---------------------------------
%     data/pw/*.mat        平面波通道数据与波束形成结果（33 个）
%     data/foc/*.mat       聚焦发射逐线结果（512 个）
%     data/noise/*.mat     噪声图像（139 个）
%     data/field/*.mat     发射声场图
%     data/res/results.mat、table_*.csv、tables.md   指标与表格
%     data/calib.mat       E0 标定（依赖阵列、脉冲与角度）
%     data/stamp_acq.mat   参数指纹
%
%   保留：data/logs（历史日志）、figures/（会重画）、data/dense 与 data/foc_rx52
%   （辅助实验与手动备份，不属于主管线），以及 data/res 里不属于本条线的文件。
%
%   为什么这些东西可以删、也必须删
%   --------------------------------
%   它们全部可以由 params + Field II 重新算出来，没有一件是"原始观测"。参数变了
%   之后，旧文件不但没用，还会成为陷阱：
%     * calib.mat 里的 dt 依赖阵列几何、角度与脉冲——留着旧 dt，图像会被整体
%       平移甚至随角度错位，而程序不会报错，只会画出不对的图；
%     * stamp_acq.mat 记录的是"这批数据由什么参数产生"，数据删了它也就失去意义，
%       留着反而会让下一次运行以为"数据还在"。
%
%   为什么 logs 与 figures 不删
%   -------------------------
%   日志是历史证据（哪次跑了多久、报了什么错），与数据是否有效无关；图会由新数据
%   重画，删不删都一样，但保留着可以和新图对照。
%
%   为什么用通配符逐类删除，而不是整个目录删掉
%   ----------------------------------------
%   data/ 下混着不属于主管线的东西：data/dense（旧的 71 角度补充实验）与
%   data/foc_rx52（修正接收孔径前的手动备份）。整目录删会顺手干掉它们，而它们
%   是复盘 B12 的证据。所以这里只按文件名模式删主管线的产物。
%
%   为什么提供 dry run
%   -----------------
%   删除不可逆。先跑 PIPELINE_RESET(P, true) 看清"要删多少、每类多少"，再决定
%   是否真的清空——尤其在本目录已经存了几百个文件的时候。
%
%   用法：run_core('cold') 会先调用本函数（非 dry run），再整条重跑。

if nargin < 2 || isempty(dryRun)
    dryRun = false;
end

patterns = { ...
    fullfile(P.dataDir, 'pw',    '*.mat'), ...
    fullfile(P.dataDir, 'foc',   '*.mat'), ...
    fullfile(P.dataDir, 'noise', '*.mat'), ...
    fullfile(P.dataDir, 'field', '*.mat'), ...
    fullfile(P.dataDir, 'res',   'results.mat'), ...
    fullfile(P.dataDir, 'res',   'table_*.csv'), ...
    fullfile(P.dataDir, 'res',   'tables.md'), ...
    fullfile(P.dataDir, 'calib.mat'), ...
    stamp_path(P)};

files = {};
for k = 1:numel(patterns)
    d = dir(patterns{k});
    for q = 1:numel(d)
        files{end+1} = fullfile(d(q).folder, d(q).name); %#ok<AGROW>
    end
end

if dryRun
    fprintf('[reset] dry run：磁盘未改动，将删除 %d 个文件\n', numel(files));
    for k = 1:numel(patterns)
        d = dir(patterns{k});
        if isempty(d), continue; end
        rel = erase(patterns{k}, [P.dataDir filesep]);
        fprintf('  %-28s %4d 个\n', rel, numel(d));
    end
    return;
end

for k = 1:numel(files)
    if exist(files{k}, 'file')
        delete(files{k});
    end
end
fprintf('[reset] 已删除 %d 个采集/分析文件，采集将从头开始\n', numel(files));
end
