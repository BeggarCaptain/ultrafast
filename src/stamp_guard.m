function stamp_guard(P, action)
%STAMP_GUARD  采集前的参数指纹闸门：参数变了就不许复用旧数据。
%
%   STAMP_GUARD(P) 检查 data/stamp_acq.mat 里记录的参数与当前 P 是否一致：
%     * 一致           —— 直接返回，采集阶段按“文件存在就跳过”继续；
%     * 没有指纹文件   —— 若已有数据则报错（旧数据来历不明），否则写入指纹；
%     * 不一致         —— 若还没有数据，重写指纹；否则报错并列出改了哪几项。
%
%   STAMP_GUARD(P, 'write') 无条件把当前参数写为指纹（供 run_core('adopt') 用）。
%
%   为什么要有这道闸门
%   ------------------
%   acquire_all_pw / acquire_all_focused / stage_noise 都是断点续跑（文件存在就
%   跳过）。没有闸门时，改了 P.elemPitch、P.angles、P.randomSeed 之类的参数再跑，程序会
%   静默复用旧数据，结果与参数对不上，而且没有任何提示。有了闸门，这种情况会在
%   采集开始前就停下，并告诉你改的是哪几项。
%
%   为什么不"警告一下接着跑"，而是直接报错
%   --------------------------------------
%   一次冷启动约 50 分钟。如果只是警告，跑到最后才发现"原来用的还是旧参数"，
%   损失的是整条时间线，而且报告里的数字会全部对不上参数表。报错把这个问题
%   提前到第 0 秒，代价只是重跑一条命令。
%
%   三种情形分别怎么处理
%   --------------------
%     一致             → 直接返回，采集按"文件存在就跳过"继续；
%     没有指纹但有数据 → 报错（数据来历不明，可能是旧版本程序留下的）；
%     不一致且没有数据 → 直接写新指纹（没有数据要保护，没必要打扰人）；
%     不一致但有数据   → 报错，并用 stamp_diff 列出改了哪几项。
%   区分"有没有数据"很重要：首次运行、或刚 pipeline_reset 之后，指纹缺失是正常
%   状态，不该报错。
%
%   版本号为什么单独查
%   ----------------
%   参数变了可以 run_core('adopt') 接受；**采集算法**变了（version +1）则不行——
%   数据格式或语义可能已经不同（例如文件里存的字段变了），旧数据根本无法解释。
%   所以 version 不一致时只给"必须 cold"一条路。
%
%   本函数只“拒绝”，不自动删数据：删除是不可逆操作，交给 run_core('cold') 或
%   pipeline_reset 显式执行。指纹文件里除了参数还记了写入时间、MATLAB 版本与
%   工程根目录，出问题时可以回溯"这份数据是什么时候、在哪台机器上产生的"。

if nargin < 2 || isempty(action)
    action = 'check';
end
f   = stamp_path(P);
cur = data_stamp(P);

% ---------------------------------------------------------- 无条件写入指纹
if strcmpi(action, 'write')
    write_stamp(f, cur, P);
    fprintf(2, '[stamp] 警告：已把当前参数登记为现有数据的来源。\n');
    fprintf(2, '        若数据其实是用旧参数算出来的，结果不可信——请改用 run_core(''cold'')。\n');
    return;
end
if ~strcmpi(action, 'check')
    error('stamp_guard:badAction', '未知 action "%s"（只支持 check / write）', action);
end

% ------------------------------------------------------------------ 比对
if ~exist(f, 'file')
    if has_data(P)
        error('stamp:missing', '%s', msg_missing(f));
    end
    write_stamp(f, cur, P);
    return;
end

old = load(f);
if ~isfield(old, 'stamp') || ~isfield(old.stamp, 'fields') || ...
        ~isfield(old.stamp, 'version')
    error('stamp:corrupt', '参数指纹文件 %s 格式不对，请删除后重跑（或用 run_core(''cold'')）。', f);
end

if old.stamp.version ~= cur.version
    error('stamp:version', ['采集算法版本变了（%d -> %d）。数据必须重算：\n' ...
        '  run_core(''cold'')  —— 约 50 min\n' ...
        '若只是改了注释或显示代码，把 src/data_stamp.m 里的 S.version 改回去即可。'], ...
        old.stamp.version, cur.version);
end

d = stamp_diff(old.stamp.fields, cur.fields);
if isempty(d)
    return;                                   % 一致：正常继续
end
if ~has_data(P)
    write_stamp(f, cur, P);                   % 没有数据要保护，直接更新指纹
    return;
end
error('stamp:mismatch', '%s', msg_mismatch(d, old));
end

% =========================================================================
function tf = has_data(P)
%HAS_DATA  是否已经存在任何采集产物。
tf = exist(fullfile(P.dataDir, 'calib.mat'), 'file') == 2 || ...
     ~isempty(dir(fullfile(P.dataDir, 'pw', '*.mat'))) || ...
     ~isempty(dir(fullfile(P.dataDir, 'foc', '*.mat'))) || ...
     ~isempty(dir(fullfile(P.dataDir, 'noise', '*.mat')));
end

function write_stamp(f, cur, P)
%WRITE_STAMP  把当前参数指纹落盘（顺带记录时间、MATLAB 版本与工程根目录）。
S = struct('stamp', cur, 'when', char(datetime('now')), ...
           'matlab', version, 'root', P.rootDir);
save(f, '-struct', 'S');
fprintf('[stamp] 参数指纹已写入 %s（版本 %d）\n', f, cur.version);
end

% =========================================================================
function s = msg_missing(f)
s = sprintf([ ...
    '[stamp] 发现采集数据，但没有参数指纹文件：\n  %s\n\n' ...
    '这些数据可能是加入指纹机制之前生成的，无法确认是否与当前参数一致。两种处理：\n' ...
    '  1) 确认数据确实由当前参数产生  ->  run_core(''adopt'')  把它登记下来\n' ...
    '  2) 不确定                      ->  run_core(''cold'')   删除旧数据并全部重算（约 50 min）\n'], f);
end

function s = msg_mismatch(d, old)
lines = strjoin(strcat({'  - '}, d), newline);
when = '';
if isfield(old, 'when'), when = old.when; end
s = sprintf([ ...
    '[stamp] 参数与现有采集数据不一致，拒绝复用。以下参数变过：\n%s\n\n' ...
    '现有数据的指纹写于 %s。处理方式（三选一）：\n' ...
    '  1) run_core(''cold'')  删除旧数据并全部重算（约 50 min）\n' ...
    '  2) run_core(''adopt'') 明确接受旧数据、只更新指纹（有风险）\n' ...
    '  3) pipeline_reset(P, true) 先看会删哪些文件，再自行处理\n'], lines, when);
end
