function acquire_all_pw(P, tag, pos, amp, keepCh)
%ACQUIRE_ALL_PW  对某个体模采集并波束形成全部平面波角度。
%
%   ACQUIRE_ALL_PW(P, TAG, POS, AMP, KEEPCH)
%
%   输入
%     P      ：参数结构体
%     tag    ：体模标记，用于文件名（'P1'、'P3' 或 'REM'）
%     pos    ：该体模的 N x 3 散射子坐标 [m]
%     amp    ：N x 1 幅度
%     keepCh ：true  ——同时保存 128 通道原始数据（E7 吞吐实验和 F6 需要）
%              false ——只保存波束形成结果（每角度省约 1.4 MB）
%
%   输出文件：data/pw/<TAG>_a###.mat，每个角度一个，内容为
%     img   ：[Nz x Nx] 该角度波束形成后的射频（single）
%     w2    ：[Nz x Nx] 插值权重平方和（single，噪声传播用）
%     t0    ：第一个采样点的绝对时刻 [s]
%     dt    ：该角度使用的时间零点校正 [s]
%     alpha ：偏转角 [rad]
%     ch    ：原始通道数据（single，仅在 keepCh 为真时保存）
%
%   循环是“断点续跑”的：文件已存在就跳过，因此中断数小时的采集后直接重跑即
%   可；每个角度的结果都在进入下一个角度之前落盘。
%
%   耗时：20 000 个散射子时约 28 s/角度（128 通道、5–60 mm），整套 11 个角度
%   约 5 分钟。

outdir = fullfile(P.dataDir, 'pw');             % 每个角度一个文件
if ~exist(outdir, 'dir'), mkdir(outdir); end

stamp_guard(P);                                 % 参数变了就拒绝复用旧数据（见 stamp_guard）

C  = load(fullfile(P.dataDir, 'calib.mat'));    % E0 标定出的 dt(alpha)
[Th, ThR] = f2_setup(P);                        % 发射句柄 + 接收句柄

tAll = tic;                                     % 总耗时，用于估算剩余时间
for i = 1:P.nAngles
    f = fullfile(outdir, sprintf('%s_a%03d.mat', tag, i));
    if exist(f, 'file'), continue; end          % 已有文件：跳过（断点续跑）

    t1 = tic;
    [rf, t0] = acquire_pw(Th, ThR, P, P.angles(i), pos, amp);  % 一次发射
    dt = C.dt(i);
    [img, w2] = das_pw(rf, t0, P.angles(i), P, P.imageX, P.imageZ, dt); % 波束形成

    S = struct('img', single(img), 'w2', single(w2), 't0', t0, ...
               'dt', dt, 'alpha', P.angles(i), 'tag', tag);
    if keepCh || (strcmp(tag,'P3') && i == P.zeroAngleIndex)   % 点目标体模的 0° 角度，留给 F6
        S.ch = single(rf);                       % 原始通道存成 single 省空间
    end
    save(f, '-struct', 'S');                     % 先落盘再做下一个

    % 打印进度与预计剩余时间（按已完成的平均速度线性外推）
    el = toc(tAll);
    fprintf('[%s] angle %2d/%d  %5.1f s  elapsed %5.1f min  eta %5.1f min\n', ...
        tag, i, P.nAngles, toc(t1), el/60, el/60*(P.nAngles-i)/i);
end
end
