function S = data_stamp(P)
%DATA_STAMP  影响“采集结果”的参数集合，用于判断缓存数据是否过期。
%
%   S = DATA_STAMP(P) 返回一个只含“会改变 data/pw、data/foc、data/noise、
%   calib.mat 内容”的那些参数的结构体。显示用、绘图用、分析用的参数
%   （displayRangeDB、figureDir、compoundCounts 等）刻意不在里面：改它们不该
%   触发重跑采集。
%
%   为什么需要它
%   ------------
%   acquire_all_pw / acquire_all_focused / stage_noise 都是“文件存在就跳过”的
%   断点续跑。好处是中断后能接着跑，坏处是：改了参数再跑，会静默复用旧数据，
%   结果与参数对不上。这个函数给出一个可比较的“参数指纹”，由 stamp_guard 在
%   每次采集前比对，不一致就报错并列出改了哪几项。
%
%   S.version
%   ---------
%   采集*算法*（而不是参数）改动时手动 +1：例如延时模型、脉冲构造、文件里存的
%   字段变了。这样即使参数没变，旧数据也会被判为过期。
%   指纹的键就是参数字段名，所以给参数改名同样算一次格式变更（v2 即此原因）。
%
%   新增参数时怎么办
%   ----------------
%   如果新参数会影响采集结果，把名字加进下面的 names 列表，并把 S.version +1；
%   如果只影响显示或分析，不要加。加错方向的后果：该重跑的没重跑（漏加）。
%
%   不跟踪的东西
%   ------------
%   * Field II / MATLAB 版本、机器硬件：指纹只看参数。换 Field II 版本或换机器后，
%     请手动 run_core('cold')，否则可能拿着旧版本算的数据继续用。
%   * 只影响显示或分析的参数（displayRangeDB、compoundCounts、figureDir 等）：
%     改了不该重跑采集。

% 采集算法版本号：改了 acquire_* / f2_setup / das_* 里任何影响数据的逻辑就 +1。
% v2：参数字段改为可读全名，指纹的键随之改变（数值一个没动）。
S.version = 2;

names = { ...
    'soundSpeed','sampleRate','centerFreq', ...      % 介质、采样率、标称频率
    'nElements','elemPitch','elemWidth','elemHeight','elemKerf', ...   % 阵列几何
    'fNumber', ...                                   % 接收 F 数（决定动态孔径）
    'pulseSigmaCycles','pulseCarrierFreq', ...       % 发射脉冲
    'txFocusDepth','txApodWidth','nTxLines','txLinePitch','nPsfLines', ...  % 聚焦发射
    'angles','nAngles', ...                          % 平面波角度
    'randomSeed','nScatterers','phantomExtentX','phantomExtentZ', ...   % 均匀体模
    'cystRadius','cystCenterX','cystCenterZ', ...    % 囊肿
    'pointTargetX','pointTargetZ','pointTargetAmp', ...  % 点目标
    'pwSnrTargetDB','nMonteCarlo', ...               % 噪声标定
    'nImageLines','pixelPitchX','pixelPitchZ','depthMin','depthMax'};   % 成像网格

S.fields = struct();
for k = 1:numel(names)
    if ~isfield(P, names{k})
        error('data_stamp:missing', ...
            '参数 P.%s 不存在。若刚删除/改名了参数，请同步更新 src/data_stamp.m 的 names 列表。', ...
            names{k});
    end
    S.fields.(names{k}) = P.(names{k});
end
end
