function stage_noise(P, kind)
%STAGE_NOISE  白色通道噪声 -> 波束形成后的噪声图像。
%
%   STAGE_NOISE(P, KIND) 为每次发射生成一份独立的噪声场，送入与信号相同的
%   DAS，保存得到的噪声图像。把噪声图像加到信号图像上就是真实的含噪帧；P1 与
%   P2（P2 = P1 - REM）共用同一份噪声实现，就像真实采集里噪声并不知道体模长
%   什么样。
%
%   KIND = 'pw'  -> data/noise/N1_a###.mat，每个平面波角度一个文件
%   KIND = 'foc' -> data/noise/NF_l###.mat，每条聚焦线一个文件
%
%   为什么噪声要加在“通道域”
%   ------------------------
%   论文的 10log10(n) SNR 增益成立的前提是：n 次发射的噪声互相独立，而信号
%   相干。若直接在波束形成后的图像上加噪声，这个前提就不成立了。这里每次发射
%   都重新取一份 randn() 噪声场，因此复合时叠加的是 n 份独立噪声（功率 x n）
%   和同一份信号（幅度 x n）。
%
%   记录长度 Nt 取自同一次发射的“信号”采集，保证噪声与信号的时间轴完全一致；
%   聚焦线并不都保存了通道数据，因此退回到同一体模的平面波记录（同阵列、同
%   深度范围，长度一致）。
%
%   断点续跑：已存在的文件直接跳过。

C = load(fullfile(P.dataDir, 'calib.mat'));
outdir = fullfile(P.dataDir, 'noise');
if ~exist(outdir, 'dir'), mkdir(outdir); end
sig = C.sigma_e;                 % stage_sigma 标定出的通道噪声标准差
Th  = f2_setup(P);               % 这里不发波，只是让 Field II 状态保持一致
rng(7, 'twister');               % 固定种子：噪声同样可复现

if strcmp(kind, 'pw')
    % ---- 每个平面波角度生成一帧噪声 ------------------------------------
    for i = 1:P.nAngles
        fs_ = fullfile(outdir, sprintf('N1_a%03d.mat', i));
        if exist(fs_, 'file'), continue; end

        S0 = load(fullfile(P.dataDir, 'pw', sprintf('P1_a%03d.mat', i)));
        Nt = size(S0.img, 1);
        if isfield(S0, 'ch'), Nt = size(S0.ch, 1); end   % 用真实的通道记录长度

        nz = sig * randn(Nt, P.N);                       % 每个阵元、每个角度独立
        [img, ~] = das_pw(nz, S0.t0, S0.alpha, P, P.x, P.z, S0.dt);

        S = struct('img', single(img), 'tag', 'N1', 'alpha', S0.alpha);
        save(fs_, '-struct', 'S');
        if mod(i,10)==0, fprintf('[noise-pw] %d/%d\n', i, P.nAngles); end
    end
else
    % ---- 每条聚焦线生成一帧噪声 ----------------------------------------
    for k = 1:P.txLines
        fs_ = fullfile(outdir, sprintf('NF_l%03d.mat', k));
        if exist(fs_, 'file'), continue; end

        S0 = load(fullfile(P.dataDir, 'foc', sprintf('P1_l%03d.mat', k)));
        % 通道记录长度：只有中心线存了 'ch'，因此退回到平面波记录
        %（同一体模、同一阵列、同一深度范围，长度必然相同）
        if isfield(S0, 'Nt')
            Nt = S0.Nt;                                  % acquire_all_focused 已存
        elseif isfield(S0, 'ch')
            Nt = size(S0.ch, 1);
        else
            Sp = load(fullfile(P.dataDir, 'pw', 'P1_a001.mat'));
            Nt = size(Sp.ch, 1);                         % 同一体模、同一距离范围
        end

        nz = sig * randn(Nt, P.N);
        [img, ~] = das_focused(nz, S0.t0, S0.xf, P, P.z, S0.dt);

        S = struct('img', single(img), 'tag', 'NF', 'xf', S0.xf);
        save(fs_, '-struct', 'S');
        if mod(k,20)==0, fprintf('[noise-foc] %d/%d\n', k, P.txLines); end
    end
end
end
