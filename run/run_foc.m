function run_foc()
%RUN_FOC  重跑聚焦采集（接收孔径修正之后）。只做采集，不做分析。
%
%   背景：改用独立接收句柄 ThR 之后，data/foc 里 10-22 的旧文件（接收通道被
%   52 阵元发射窗加权）作废，需要重新采集 P3 / P1 / REM 的全部扫描线。
%   平面波不受影响——收发变迹都是单位窗，实测新旧数据逐位相同，所以 data/pw
%   保持不动。
%
%   用法（后台跑，不占用交互会话）：
%       J = batch(@run_foc, 0, {}, 'Profile', 'Processes');
%       diary(J)            % 看进度
%
%   耗时参考：P1 128 线约 15 min，REM 约 1 min，P3 256 线约 1 min。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);

fprintf('=== focused re-acquire start %s ===\n', char(datetime('now')));

% P3：8 个点目标，PSF 用 256 条细线
[pp, pa] = phantom_points(P);
acquire_all_focused(P, 'P3', pp, pa, P.psfLineX, false);

% P1：均匀散斑体模，保留中心线的原始通道（F6 与 E7 要用）
[pos, amp] = phantom_uniform(P);
acquire_all_focused(P, 'P1', pos, amp, P.txLineX, true);

% REM：囊肿内部的散射子，用于合成无回声体模 P2 = P1 - REM
[rp, ra] = phantom_cyst_removed(P, pos, amp);
acquire_all_focused(P, 'REM', rp, ra, P.txLineX, false);

fprintf('=== focused re-acquire done %s ===\n', char(datetime('now')));
end
