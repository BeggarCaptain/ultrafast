function stage_sigma(P)
%STAGE_SIGMA  标定通道噪声水平，并用蒙特卡洛验证解析噪声传播公式（设计 4.5 节）。
%
%   目标：单个平面波在 z = 30 mm 处的 SNR 必须等于 P.snrPW30（20 dB）。之后
%   所有方法、所有体模共用同一个 sigma_e，SNR 比较才是公平的。
%
%   这里做两件事：
%     1. 由实测信号功率与实测波束形成权重和反解 sigma_e：
%        sigma_e = sqrt(Ps / (10^(SNR/10) * W))。
%     2. 在小 ROI 上把解析噪声增益 sigma_N^2 = sigma_e^2 * sum(w^2) 与蒙特卡洛
%        结果对比。两者若不一致，报告里整节 SNR 结论都作废。
%
%   输出：sigma_e 与蒙特卡洛比值追加写入 calib.mat。

fprintf('\n===== noise calibration =====\n');
C = load(fullfile(P.dataDir,'calib.mat'));
S = load(fullfile(P.dataDir,'pw',sprintf('P1_a%03d.mat', P.iZero)));   % 0° 平面波、P1 体模
img = double(S.img);  w2 = double(S.w2);            % 波束形成信号与权重和

% ---- ROI：30 mm 附近 10 x 10 mm，即发射焦点附近 ------------------------
[XX, ZZ] = meshgrid(P.x, P.z);
roi = abs(XX) <= 5e-3 & ZZ >= 25e-3 & ZZ <= 35e-3;

Ps = mean(img(roi).^2);                             % ROI 内平均信号功率
W  = mean(w2(roi));                                 % 平均权重平方和
% 要求 10*log10(Ps / (sigma_e^2 W)) = snrPW30，反解 sigma_e
sigma_e = sqrt(Ps / (10^(P.snrPW30/10) * W));
fprintf('signal power %.3e, mean sum(w^2) %.1f -> sigma_e %.3e\n', Ps, W, sigma_e);

% ---------------------------------------------------- 蒙特卡洛验证
% 用小网格是为了让 20 次重复足够便宜。每次生成独立的白色噪声场，送入同一个
% DAS，再把经验标准差与 sigma_e*sqrt(w2) 比较。
xr = (-5:0.1711:5)*1e-3;  zr = (25:0.1:35)*1e-3;
Nt = size(S.ch, 1);
[~, w2r] = das_pw(zeros(Nt,P.N,'single'), S.t0, 0, P, xr, zr, S.dt);  % 只要权重
acc = zeros(numel(zr), numel(xr));  acc2 = acc;
rng(11,'twister');                                  % 固定种子，验证可复现
for r = 1:P.nMC
    nz = sigma_e * randn(Nt, P.N);                  % 白色通道噪声
    nr = das_pw(nz, S.t0, 0, P, xr, zr, S.dt);      % 波束形成后的噪声
    acc  = acc  + nr;
    acc2 = acc2 + nr.^2;
end
emp = sqrt(max(acc2/P.nMC - (acc/P.nMC).^2, 0));    % 经验标准差（E[x^2]-E[x]^2）
ana = sigma_e * sqrt(w2r);                          % 解析标准差
ratio = emp ./ ana;
fprintf('Monte-Carlo (%d reps): empirical/analytic noise ratio = %.3f +- %.3f\n', ...
    P.nMC, mean(ratio(:)), std(ratio(:)));

% ---- 追加写回 calib.mat（其余字段保持不变） ----------------------------
C.sigma_e = sigma_e;
C.mcRatio = mean(ratio(:));
C.mcRatioStd = std(ratio(:));
C.Ps30 = Ps;  C.W30 = W;
save(fullfile(P.dataDir,'calib.mat'), '-struct', 'C');
fprintf('sigma_e saved to calib.mat\n');
end
