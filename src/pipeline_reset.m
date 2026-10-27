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
