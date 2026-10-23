function run_dense_snr()
%RUN_DENSE_SNR  补充实验：把散射子密度提高 4 倍，能否拿回理想的 10log10(n) 增益？
%
%   动机
%   ----
%   主实验测得 71 角度的 SNR 增益是 15.6 dB，而理想值 18.5 dB。第一个怀疑是
%   二维体模太稀（每个分辨单元只有约 2 个散射子），散斑没完全发展。本脚本直接
%   检验这个假设。
%
%   设计
%   ----
%   同样的阵列、同样的脉冲、同样的 sigma_e、同样的 71 个角度，但把两万个散射子
%   压缩到 20 x 30 mm 内（33 个/mm²，约为 P1 的 4 倍）。只需要平面波发射（不需要
%   聚焦参照），因此约 30 分钟即可跑完。
%
%   输出
%   ----
%   data/dense/D_a###.mat        每个角度的信号与噪声图像
%   data/res/dense_snr.mat       各档位在 30 mm 处的相对 SNR
%
%   结论（见报告 5.9）：增益只提高 0.04–0.37 dB，所以密度不是缺额的原因；真正
%   的原因在散斑角度间的相干性（报告 6.5）。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'run'));
P = params(root);
C = load(fullfile(P.dataDir, 'calib.mat'));

outdir = fullfile(P.dataDir, 'dense');
if ~exist(outdir, 'dir'), mkdir(outdir); end
logdir = fullfile(P.dataDir, 'logs');
diary(fullfile(logdir, 'dense.log')); diary on;
fprintf('=== dense phantom SNR check %s ===\n', char(datetime('now')));
t0 = tic;

try
    % ---- 致密体模：散射子个数不变，密度提高 4 倍 -----------------------
    % 刻意用与 P1 不同的种子：这是一次独立实现，而不是主体模的子集。
    rng(P.seed + 1, 'twister');
    n = P.scatN;
    pos = [ -10e-3 + 20e-3*rand(n,1), zeros(n,1), 15e-3 + 30e-3*rand(n,1) ];
    amp = randn(n,1);
    fprintf('density %.1f /mm2\n', n/(20*30));

    % ---- 逐角度采集，每次配一份独立噪声 ------------------------------
    % 信号与噪声都在这里波束形成，这样 SNR 的计算方式与主分析完全一致。
    Th = f2_setup(P);
    for i = 1:P.nAngles
        f = fullfile(outdir, sprintf('D_a%03d.mat', i));
        if exist(f, 'file'), continue; end            % 断点续跑
        [rf, t00] = acquire_pw(Th, P, P.angles(i), pos, amp);
        [img, w2] = das_pw(rf, t00, P.angles(i), P, P.x, P.z, C.dt(i));
        nz  = C.sigma_e * randn(size(rf));            % 与主实验相同的 sigma_e
        niv = das_pw(nz, t00, P.angles(i), P, P.x, P.z, C.dt(i));
        S = struct('img', single(img), 'nz', single(niv), 'w2', single(w2), ...
                   't0', t00, 'dt', C.dt(i));
        save(f, '-struct', 'S');
        if mod(i,10) == 0
            fprintf('  angle %2d/%d  elapsed %.1f min\n', i, P.nAngles, toc(t0)/60);
        end
    end

    % ---- 各档位的 SNR 增益 --------------------------------------------
    % 定义与主分析相同：5x5 mm 窗口内的 10log10(<s^2>/<n^2>)，参考点为 30 mm 处
    % 的单角度值。
    sets = subset_indices(P);
    zdep = (15:5:45)*1e-3;
    gain = zeros(1, numel(P.subsetN));
    [zz, s0] = snr_curve_(load(fullfile(outdir,sprintf('D_a%03d.mat', P.iZero))), P);
    ref = interp1(zz, s0, 30e-3);                     % 0 dB 参考
    for q = 1:numel(P.subsetN)
        S = zeros(numel(P.z), numel(P.x), 'single');
        N = zeros(numel(P.z), numel(P.x), 'single');
        for i = sets{q}
            T = load(fullfile(outdir, sprintf('D_a%03d.mat', i)));
            S = S + T.img;  N = N + T.nz;             % 相干求和，与生产代码一致
        end
        [zz, s] = snr_curve_rm(S, N, P);
        gain(q) = interp1(zz, s, 30e-3) - ref;
        fprintf('  n=%2d: relative SNR at 30 mm = %5.2f dB (theory %.2f)\n', ...
            P.subsetN(q), gain(q), 10*log10(P.subsetN(q)));
    end
    save(fullfile(P.dataDir,'res','dense_snr.mat'), 'P', 'gain', 'zdep');
    fprintf('=== dense check done in %.1f min ===\n', toc(t0)/60);
catch err
    fprintf(2, 'ERROR: %s\n', err.message);
end
diary off;
end

function [zz, snr] = snr_curve_(S, P)
%SNR_CURVE_  包装函数：对一份已波束形成的角度文件求 SNR 曲线。
[zz, snr] = snr_curve_rm(S.img, S.nz, P);
end

function [zz, snr] = snr_curve_rm(rf, nrf, P)
%SNR_CURVE_RM  10log10(信号功率/噪声功率) 的深度曲线。
%   两者都用可分离的 5x5 mm 盒平滑；沿 x 取中位数，每个深度得到一个稳健值。
ps = double(rf).^2;  pn = double(nrf).^2;

h1 = ones(round(5e-3/P.dz),1); h1 = h1/sum(h1);       % 轴向 5 mm 盒
h2 = ones(1, round(5e-3/P.dx));  h2 = h2/sum(h2);     % 横向 5 mm 盒
ps = imfilter(imfilter(ps, h1, 'replicate'), h2, 'replicate');
pn = imfilter(imfilter(pn, h1, 'replicate'), h2, 'replicate');

S  = 10*log10(ps ./ max(pn, realmin));
zz = P.z;
snr = median(S, 2);
end
