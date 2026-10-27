function run_E7_gpu_throughput()
%RUN_E7_GPU_THROUGHPUT  E7（可选）：算法处理帧率。
%
%   这里测什么
%   ----------
%   把“一帧 11 角度的通道数据”处理成一张 B 模式图，在 CPU 上和在 GPU 上各需要
%   多少墙钟时间，并换算成“算法处理帧率”（每秒能处理多少帧）。
%
%   这不是什么
%   ----------
%   它不是采集帧率 f_PRF/n。采集帧率由发射序列决定，和机器快慢毫无关系
%   （见设计文档 7.1 节）。两者回答的是不同问题，报告里也分开列两张表。
%
%   Field II 在这里完全不参与：通道数据就是核心实验已经存下的
%   data/pw/P1_a###.mat（keepCh = true）。这样 CPU 与 GPU 版本拿到的输入逐字节
%   完全相同。
%
%   计时协议
%   --------
%   * 先跑 nWarm 次预热并丢弃（内存分配、JIT、GPU 升频都在这一阶段）；
%   * 正式重复 nRep 次，用 tic/toc 计时，取中位数；
%   * GPU 给两个数：
%       compute-only ：数据已在显存里（纯算法速度），
%       end-to-end   ：含一帧的 H2D 传输（把 PCIe 带宽也算进去）；
%   * 只有当 GPU 与 CPU 图像的 NRMSE 小于 1e-3 时，才认为速度数字可信。
%
%   运行条件：CUDA GPU + Parallel Computing Toolbox。没有的话脚本打印
%   “未执行”后返回，项目其他部分都不依赖它。
%
%   输出：data/res/e7_throughput.mat、figures/F15_e7_throughput.png
%
%   注意：本机没有 CUDA 设备（gpuDeviceCount = 0），因此报告中的 E7 未执行；
%   脚本已就绪，换到有 GPU 的机器上可直接运行。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);
C = load(fullfile(P.dataDir, 'calib.mat'));

outdir = fullfile(P.dataDir, 'res');
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ---------------------------------------------------------------- GPU 检查
% 明确记录设备型号：没有设备名的速度数字是不可复现的。
hasGPU = false;  gpuInfo = 'none';
try
    if gpuDeviceCount > 0
        g = gpuDevice;
        hasGPU = true;
        gpuInfo = sprintf('%s (%.1f GB)', g.Name, g.TotalMemory/1e9);
    end
catch
    % 没有 GPU，或 Parallel Computing Toolbox 未授权：继续往下走
end
fprintf('GPU: %s\n', gpuInfo);

% -------------------------------------------------------------- 读入一帧
% 11 个通道记录（每个约 1.2 MB，single）构成一帧的输入。
chs = cell(1, P.nAngles);
t0s = zeros(1, P.nAngles);  dts = zeros(1, P.nAngles);
for i = 1:P.nAngles
    S = load(fullfile(P.dataDir, 'pw', sprintf('P1_a%03d.mat', i)));
    chs{i} = S.ch;  t0s(i) = S.t0;  dts(i) = S.dt;
end
fprintf('frame size: %d x %d per angle, %d angles\n', size(chs{1}), P.nAngles);

nList = P.gpuBench.counts;  nWarm = P.gpuBench.warmup;  nRep = P.gpuBench.repeat;
res = struct('gpu', gpuInfo, 'nList', nList, ...
             'cpu', nan(1,numel(nList)), 'gpuC', nan(1,numel(nList)), ...
             'gpuE', nan(1,numel(nList)), 'nrmse', nan(1,numel(nList)));

% ------------------------------------------------------------------- CPU
% 基线：生产代码路径，单线程，不加项目里本来没有的向量化技巧。
for q = 1:numel(nList)
    idx = subset_idx(P, nList(q));
    f = @() pipeline_cpu(chs, idx, t0s, dts, P);
    for w = 1:nWarm, f(); end                 % 预热，丢弃
    tt = zeros(1, nRep);
    for r = 1:nRep, tic; f(); tt(r) = toc; end
    res.cpu(q) = median(tt);                  % 中位数：对系统抖动更稳
    fprintf('CPU n=%2d: %6.1f ms  (%5.1f fps)\n', nList(q), res.cpu(q)*1e3, 1/res.cpu(q));
    if nList(q) == P.nAngles, imgCPU = f(); end   % 作为 NRMSE 校验的参考图
end

% ------------------------------------------------------------------- GPU
if hasGPU
    for q = 1:numel(nList)
        idx = subset_idx(P, nList(q));

        % ---- compute-only：这一帧已经在显存里 --------------------------
        gch = cellfun(@(x) gpuArray(x), chs(idx), 'UniformOutput', false);
        f = @() pipeline_gpu(gch, idx, t0s, dts, P);
        for w = 1:nWarm, f(); end
        tt = zeros(1, nRep);
        for r = 1:nRep, tic; f(); tt(r) = toc; end
        res.gpuC(q) = median(tt);

        % ---- end-to-end：把主机到显存的传输算进去 ----------------------
        f2 = @() pipeline_gpu_e2e(chs, idx, t0s, dts, P);
        for w = 1:nWarm, f2(); end
        tt = zeros(1, nRep);
        for r = 1:nRep, tic; f2(); tt(r) = toc; end
        res.gpuE(q) = median(tt);

        % ---- 对完整帧（全部角度）做一致性校验 ---------------------------
        if nList(q) == P.nAngles
            imgGPU = gather(pipeline_gpu(gch, idx, t0s, dts, P));
            res.nrmse(q) = norm(double(imgGPU(:)) - double(imgCPU(:))) / norm(double(imgCPU(:)));
        end
        fprintf('GPU n=%2d: compute %6.1f ms (%5.1f fps), end-to-end %6.1f ms\n', ...
            nList(q), res.gpuC(q)*1e3, 1/res.gpuC(q), res.gpuE(q)*1e3);
    end
    res.speedup = res.cpu ./ res.gpuC;
    % 跑得快但结果不对的 GPU 数字，绝不能进报告。
    assert(res.nrmse(end) < 1e-3, 'GPU/CPU mismatch, speed numbers not reportable');
else
    fprintf('E7 not executed: no CUDA GPU available on this machine.\n');
end

save(fullfile(outdir, 'e7_throughput.mat'), '-struct', 'res');
if hasGPU
    plot_e7(res, fullfile(P.figureDir, 'F15_e7_throughput.png'));
end
end

% =========================================================================
function img = pipeline_cpu(chs, idx, t0s, dts, P)
%PIPELINE_CPU  CPU 上的一帧：DAS -> 相干求和 -> 包络 -> dB。
%   这就是生产算法本身，因此计时有代表性，而不是某个内核的微基准。
S = zeros(numel(P.imageZ), numel(P.imageX), 'single');
for j = 1:numel(idx)
    i = idx(j);
    S = S + das_pw(chs{i}, t0s(i), P.angles(i), P, P.imageX, P.imageZ, dts(i));
end
A = abs(hilbert(double(S)));                       % 求和之后再取包络
img = 20*log10(max(A, realmin) / max(A(:)));       % 对数压缩，峰值为 0 dB
end

function img = pipeline_gpu(gch, idx, t0s, dts, P)
%PIPELINE_GPU  同样的流水线，全部数组都在显存里。
%   hilbert 直接在 gpuArray 上执行，包络留在显存，直到调用者 gather 才过总线。
S = gpuArray.zeros(numel(P.imageZ), numel(P.imageX), 'single');
for j = 1:numel(idx)
    i = idx(j);
    S = S + das_pw_gpu(gch{j}, t0s(i), P.angles(i), P, dts(i));
end
A = abs(hilbert(S));
img = 20*log10(max(A, realmin) / max(A(:)));
end

function img = pipeline_gpu_e2e(chs, idx, t0s, dts, P)
%PIPELINE_GPU_E2E  端到端版本：H2D 传输也计入计时。
%   列出来是为了把“算法速度”和“PCIe 带宽”分开：一帧约 85 MB，传输不可忽略。
gch = cellfun(@(x) gpuArray(x), chs(idx), 'UniformOutput', false);
img = pipeline_gpu(gch, idx, t0s, dts, P);
end

function idx = subset_idx(P, n)
%SUBSET_IDX  角度数为 n 的那个复合档位的下标集合。
S = subset_indices(P);
idx = S{P.compoundCounts == n};
end

function plot_e7(res, fname)
%PLOT_E7  图 F15：单帧耗时与算法处理帧率随 n 的变化。
%   右图的标题特意写明“不是采集帧率”，避免报告里把两者混为一谈。
fig = figure('Position', [100 100 900 380]);
subplot(1,2,1);
semilogy(res.nList, res.cpu*1e3, 'o-', 'LineWidth', 1.3); hold on;
semilogy(res.nList, res.gpuC*1e3, 's-', 'LineWidth', 1.3);
semilogy(res.nList, res.gpuE*1e3, 'd--', 'LineWidth', 1.1);
xlabel('n angles'); ylabel('time per frame (ms)'); grid on;
legend('CPU','GPU compute-only','GPU end-to-end','Location','northwest');
title('Algorithmic processing time');
subplot(1,2,2);
plot(res.nList, 1./res.cpu, 'o-', 'LineWidth', 1.3); hold on;
plot(res.nList, 1./res.gpuC, 's-', 'LineWidth', 1.3);
xlabel('n angles'); ylabel('algorithmic frame rate (fps)'); grid on;
legend('CPU','GPU','Location','northeast');
title('NOT the acquisition frame rate');
save_fig(fig, fname); close(fig);
end
