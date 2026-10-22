function run_E7_gpu_throughput()
%RUN_E7_GPU_THROUGHPUT  E7 (optional): algorithmic processing frame rate.
%
%   Measures the wall-clock time needed to turn one frame of 71-angle channel
%   data into a B-mode image, on the CPU and (if available) on a GPU, and
%   converts it into an *algorithmic* frame rate.  This is NOT the acquisition
%   frame rate f_PRF/n (design 7.1) and Field II is never re-run here: the
%   channel data are the ones stored by the core run.
%
%   Requirements: CUDA GPU + Parallel Computing Toolbox.  Without them the
%   script records "not executed" and returns without touching the report.
%
%   Output: data/res/e7_throughput.mat and figures/F15_e7_throughput.png

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);
C = load(fullfile(P.dataDir, 'calib.mat'));

outdir = fullfile(P.dataDir, 'res');
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ---------------------------------------------------------------- GPU check
hasGPU = false;  gpuInfo = 'none';
try
    if gpuDeviceCount > 0
        g = gpuDevice;
        hasGPU = true;
        gpuInfo = sprintf('%s (%.1f GB)', g.Name, g.TotalMemory/1e9);
    end
catch
end
fprintf('GPU: %s\n', gpuInfo);

% -------------------------------------------------------------- load a frame
chs = cell(1, P.nAngles);
t0s = zeros(1, P.nAngles);  dts = zeros(1, P.nAngles);
for i = 1:P.nAngles
    S = load(fullfile(P.dataDir, 'pw', sprintf('P1_a%03d.mat', i)));
    chs{i} = S.ch;  t0s(i) = S.t0;  dts(i) = S.dt;
end
fprintf('frame size: %d x %d per angle, %d angles\n', size(chs{1}), P.nAngles);

nList = P.e7.nList;  nWarm = P.e7.warmup;  nRep = P.e7.repeat;
res = struct('gpu', gpuInfo, 'nList', nList, ...
             'cpu', nan(1,numel(nList)), 'gpuC', nan(1,numel(nList)), ...
             'gpuE', nan(1,numel(nList)), 'nrmse', nan(1,numel(nList)));

% ------------------------------------------------------------------- CPU
for q = 1:numel(nList)
    idx = subset_idx(P, nList(q));
    f = @() pipeline_cpu(chs, idx, t0s, dts, P);
    for w = 1:nWarm, f(); end
    tt = zeros(1, nRep);
    for r = 1:nRep, tic; f(); tt(r) = toc; end
    res.cpu(q) = median(tt);
    fprintf('CPU n=%2d: %6.1f ms  (%5.1f fps)\n', nList(q), res.cpu(q)*1e3, 1/res.cpu(q));
    if nList(q) == 71, imgCPU = f(); end
end

% ------------------------------------------------------------------- GPU
if hasGPU
    for q = 1:numel(nList)
        idx = subset_idx(P, nList(q));
        gch = cellfun(@(x) gpuArray(x), chs(idx), 'UniformOutput', false);
        f = @() pipeline_gpu(gch, idx, t0s, dts, P);
        for w = 1:nWarm, f(); end
        tt = zeros(1, nRep);
        for r = 1:nRep, tic; f(); tt(r) = toc; end
        res.gpuC(q) = median(tt);

        f2 = @() pipeline_gpu_e2e(chs, idx, t0s, dts, P);
        for w = 1:nWarm, f2(); end
        tt = zeros(1, nRep);
        for r = 1:nRep, tic; f2(); tt(r) = toc; end
        res.gpuE(q) = median(tt);

        if nList(q) == 71
            imgGPU = gather(pipeline_gpu(gch, idx, t0s, dts, P));
            res.nrmse(q) = norm(double(imgGPU(:)) - double(imgCPU(:))) / norm(double(imgCPU(:)));
        end
        fprintf('GPU n=%2d: compute %6.1f ms (%5.1f fps), end-to-end %6.1f ms\n', ...
            nList(q), res.gpuC(q)*1e3, 1/res.gpuC(q), res.gpuE(q)*1e3);
    end
    res.speedup = res.cpu ./ res.gpuC;
    assert(res.nrmse(5) < 1e-3, 'GPU/CPU mismatch, speed numbers not reportable');
else
    fprintf('E7 not executed: no CUDA GPU available on this machine.\n');
end

save(fullfile(outdir, 'e7_throughput.mat'), '-struct', 'res');
if hasGPU
    plot_e7(res, fullfile(P.figDir, 'F15_e7_throughput.png'));
end
end

% =========================================================================
function img = pipeline_cpu(chs, idx, t0s, dts, P)
S = zeros(numel(P.z), numel(P.x), 'single');
for j = 1:numel(idx)
    i = idx(j);
    S = S + das_pw(chs{i}, t0s(i), P.angles(i), P, P.x, P.z, dts(i));
end
A = abs(hilbert(double(S)));
img = 20*log10(max(A, realmin) / max(A(:)));
end

function img = pipeline_gpu(gch, idx, t0s, dts, P)
S = gpuArray.zeros(numel(P.z), numel(P.x), 'single');
for j = 1:numel(idx)
    i = idx(j);
    S = S + das_pw_gpu(gch{j}, t0s(i), P.angles(i), P, dts(i));
end
A = abs(hilbert(S));
img = 20*log10(max(A, realmin) / max(A(:)));
end

function img = pipeline_gpu_e2e(chs, idx, t0s, dts, P)
gch = cellfun(@(x) gpuArray(x), chs(idx), 'UniformOutput', false);
img = pipeline_gpu(gch, idx, t0s, dts, P);
end

function idx = subset_idx(P, n)
S = subset_indices(P);
idx = S{P.subsetN == n};
end

function plot_e7(res, fname)
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
