function run_test_p3()
%RUN_TEST_P3  点目标体模（P3）上的端到端快速自检。
%
%   这是进入数小时 P1 采集之前必须通过的冒烟测试：它会跑 E0（若 calib.mat
%   不存在）、用两种发射模式采集 8 个点目标、测量 PSF，并输出一张两联的调试图。
%   总耗时约 5 分钟。
%
%   什么样的结果算“通过”
%   --------------------
%     * 两种模式都把 30 mm 处的目标重建在 30.00 ± 0.05 mm（否则时间零点或延时
%       模型有问题）；
%     * 单平面波 -10 dB 横向约 1.0 mm；聚焦约 0.8 mm；轴向约 0.18 mm；
%     * 聚焦线边缘/中心通道能量比很小，说明 52 阵元发射变迹确实生效了。
%
%   Field II 时间轴的问题（DAS 里漏掉 +x1 sin(alpha)/c）在这里表现为几个毫米的
%   横向偏移；在 8 个散射子上排查，比在两万个散射子上排查容易得多。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);

% ---- 同时写日志，便于日后追溯 ------------------------------------------
if ~exist(fullfile(P.dataDir,'logs'),'dir'), mkdir(fullfile(P.dataDir,'logs')); end
diary(fullfile(P.dataDir,'logs','test_p3.log')); diary on;
fprintf('=== P3 test %s ===\n', datestr(now));

% ---- E0 是前提：DAS 需要 dt --------------------------------------------
if ~exist(fullfile(P.dataDir,'calib.mat'),'file')
    stage0_calib(P);
end

% ---- 采集点目标（71 个平面波 + 256 条聚焦线） --------------------------
% 只有 8 个散射子，两种采集都远不到一分钟。
[pp, pa] = phantom_points(P);
acquire_all_pw(P, 'P3', pp, pa, false);
acquire_all_focused(P, 'P3', pp, pa, P.xpsf, false);

% ---- 拼接与测量 --------------------------------------------------------
Sp = load(fullfile(P.dataDir,'pw',sprintf('P3_a%03d.mat', P.iZero)));   % 中心角度（alpha = 0）
env = abs(hilbert(double(Sp.img)));
m = metrics_resolution(env, P.x, P.z, 0, 30e-3);
fprintf('PW 0 deg : peak (%.2f, %.2f) mm, w10 %.2f mm, w6 %.2f mm\n', ...
    m.xpk*1e3, m.zpk*1e3, m.w10*1e3, m.w6*1e3);

foc = zeros(numel(P.z), numel(P.xpsf));              % 拼 256 条线
for k = 1:numel(P.xpsf)
    S = load(fullfile(P.dataDir,'foc',sprintf('P3_l%03d.mat',k)));
    foc(:,k) = S.img;
end
envF = abs(hilbert(double(foc)));
mF = metrics_resolution(envF, P.xpsf, P.z, 0, 30e-3);
fprintf('focused  : peak (%.2f, %.2f) mm, w10 %.2f mm, w6 %.2f mm\n', ...
    mF.xpk*1e3, mF.zpk*1e3, mF.w10*1e3, mF.w6*1e3);

% ---- 变迹检查 ----------------------------------------------------------
% 聚焦发射只用 52 个阵元。如果变迹没生效，边缘阵元的能量会和中心一样大，等效
% 发射 F 数就会远小于 1.75。
Sf = load(fullfile(P.dataDir,'foc','P3_l128.mat'));  % 中心线，存了 'ch'
if isfield(Sf,'ch')
    e = mean(abs(Sf.ch),1);
    fprintf('transmit apodisation: edge/centre channel energy ratio %.3f\n', ...
        mean(e([1:3 end-2:end]))/mean(e(60:68)));
end

% ---- 调试图：两种模式并排 ----------------------------------------------
% 文件名带下划线前缀，明确它是诊断图，不属于报告插图 F1–F13。
fig = figure('Position',[100 100 900 420]);
subplot(1,2,1);
imagesc(P.x*1e3, P.z*1e3, to_db(env,P.drDB)); axis image; colormap(gca,gray);
caxis([-P.drDB 0]); xlabel('x (mm)'); ylabel('z (mm)'); title('PW, 0 deg');
subplot(1,2,2);
imagesc(P.xpsf*1e3, P.z*1e3, to_db(envF,P.drDB)); axis image; colormap(gca,gray);
caxis([-P.drDB 0]); xlabel('x (mm)'); ylabel('z (mm)'); title('focused, 30 mm');
save_fig(fig, fullfile(P.figDir,'_debug_p3.png'), 150);

fprintf('=== P3 test done ===\n');
diary off;
end
