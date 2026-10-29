function acquire_all_focused(P, tag, pos, amp, xf, keepCh)
%ACQUIRE_ALL_FOCUSED  对某个体模采集并波束形成全部聚焦扫描线。
%
%   ACQUIRE_ALL_FOCUSED(P, TAG, POS, AMP, XF, KEEPCH)
%
%   输入
%     P      ：参数结构体
%     tag    ：体模标记，用于文件名（'P1'、'P3' 或 'REM'）
%     pos    ：N x 3 散射子坐标 [m]
%     amp    ：N x 1 幅度
%     xf     ：线的横向位置向量 [m]——散斑/囊肿体模用 P.txLineX（128 线、0.33 mm），
%              点目标体模用 P.psfLineX（256 线、0.171 mm）
%     keepCh ：true 时额外保存中心线那一线的原始通道（否则每线多 4.5 MB，而
%              画图只需要一条线）
%
%   输出文件：data/foc/<TAG>_l###.mat，每条线一个，内容为
%     img ：[Nz x 1] 该线波束形成后的射频（single）
%     w2  ：[Nz x 1] 插值权重平方和
%     t0  ：第一个采样点的绝对时刻 [s]
%     dt  ：全局时间零点校正（E0 得到的 dt_foc）
%     xf  ：线位置 [m]
%     Nt  ：通道记录长度（噪声阶段需要）
%
%   每条线发射一次（发射端不做合成孔径），这就是作为参照方法的常规聚焦 B 模式。
%   为什么不在这里做"发射合成孔径"：那会变成另一种方法（合成孔径聚焦），就不再
%   是论文要对照的"常规聚焦 B 模式"了。它的价值正在于代表另一端——每帧 128 次
%   发射、焦区质量最好、焦外退化。
%   断点续跑方式与 acquire_all_pw 相同。
%
%   耗时：20 000 个散射子时 17–28 s/线（随机器负载），128 线约 58 min。早期版本
%   是 7 s/线：那时收发共用一个句柄，Field II 能复用空间冲激响应；改用独立接收
%   句柄后无法复用，单线变慢——这是修正接收孔径（B12）的必要代价。

outdir = fullfile(P.dataDir, 'foc');
if ~exist(outdir, 'dir'), mkdir(outdir); end

stamp_guard(P);                                 % 参数变了就拒绝复用旧数据（见 stamp_guard）

C  = load(fullfile(P.dataDir, 'calib.mat'));    % E0 得到的 dt_foc
[Th, ThR] = f2_setup(P);                        % 发射句柄 + 接收句柄

nL   = numel(xf);
tAll = tic;
for k = 1:nL
    f = fullfile(outdir, sprintf('%s_l%03d.mat', tag, k));
    if exist(f, 'file'), continue; end          % 已有文件：跳过（断点续跑）

    t1 = tic;
    [rf, t0] = acquire_focused(Th, ThR, P, xf(k), pos, amp);  % 一次聚焦发射
    [img, w2] = das_focused(rf, t0, xf(k), P, P.imageZ, C.dt_foc); % 只重建这一条线

    S = struct('img', single(img), 'w2', single(w2), 't0', t0, ...
               'dt', C.dt_foc, 'xf', xf(k), 'tag', tag, 'Nt', size(rf,1));
    if keepCh && k == round(nL/2)
        S.ch = single(rf);                      % 只存中心线，供画图用
    end
    save(f, '-struct', 'S');

    % 每 10 条线打印一次：128 条线逐条打印会刷屏
    el = toc(tAll);
    if mod(k, 10) == 0 || k == nL
        fprintf('[%s] line %3d/%d  %5.1f s  elapsed %5.1f min  eta %5.1f min\n', ...
            tag, k, nL, toc(t1), el/60, el/60*(nL-k)/k);
    end
end
end
