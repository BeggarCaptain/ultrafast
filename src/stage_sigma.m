function stage_sigma(P)
%STAGE_SIGMA  标定通道噪声水平，并用蒙特卡洛验证解析噪声传播公式（设计 4.5 节）。
%
%   它在整条链路里的位置
%   --------------------
%   采集（P1 的 0° 平面波）→ 本函数给出一个数 sigma_e（每个通道、每个采样点的
%   噪声标准差）→ stage_noise 用它生成噪声通道、波束形成出噪声图像；同时，所有
%   方法、所有深度的 SNR 都可以用解析式 sigma_N^2 = sigma_e^2 * sum(w^2) 直接算，
%   不必每次真跑一遍噪声。
%
%   为什么要在"图像上"定一个 SNR 目标
%   ---------------------------------
%   通道级的信噪比无法比较不同方法：聚焦发射与平面波发射的发射能量不同、接收
%   孔径随深度变化，通道上的数会跟着变。真正该固定的是**图像上**的信噪比，所以
%   这里规定：单个平面波在 z = 30 mm（发射焦点、也是复合公共区最宽处）的 SNR
%   等于 P.pwSnrTargetDB = 20 dB。定了这一个点，sigma_e 就唯一确定，其余方法、
%   其余深度的 SNR 都是结果而不是输入。
%
%   怎么反解 sigma_e
%   ----------------
%   DAS 对每个像素是一次线性加权求和（权重含孔径掩膜与插值），若输入噪声是
%   方差 sigma_e^2 的独立同分布白噪声，输出噪声方差就是 sigma_e^2 * w2，其中
%   w2 = sum(w0^2 + w1^2) 由 das_pw 一并算出。于是图像 SNR 为
%       SNR = Ps / (sigma_e^2 * W)          Ps = 该 ROI 的信号功率，W = 平均 w2
%   要求它等于 10^(20/10)，反解
%       sigma_e = sqrt(Ps / (10^(SNR/10) * W))
%   注意 Ps 是**波束形成之后**的功率、W 是同一像素上的权重和：两者相除正好抵掉
%   DAS 的相干增益，所以得到的是"每个通道"的噪声，换一种方法（例如 11 角度复合）
%   可以直接套用同一个 sigma_e。
%
%   为什么还要蒙特卡洛验证
%   ----------------------
%   上面的推导依赖三个假设：通道噪声独立同分布、DAS 的插值权重是线性的、孔径
%   掩膜在像素间不相关。任何一个不成立，解析式就会系统性偏大或偏小，而报告里
%   整节 SNR 结论都建立在它之上。所以用 20 次独立的噪声实现跑真实 DAS，把经验
%   标准差与解析值逐像素相比；比值应接近 1（本机 0.967 ± 0.157），偏差在 3.3 %
%   以内，说明解析式可用。
%
%   为什么 ROI 取 30 mm 附近 10 x 10 mm
%   -----------------------------------
%   那是发射焦点、信号最强、散斑统计最平稳的一段；ROI 太小（例如 2 x 2 mm）样本
%   不足，w2 与信号功率的估计都会被几个像素主导。
%
%   输出：sigma_e、蒙特卡洛比值追加写入 calib.mat（其余字段保持不变）。

fprintf('\n===== noise calibration =====\n');
C = load(fullfile(P.dataDir,'calib.mat'));
S = load(fullfile(P.dataDir,'pw',sprintf('P1_a%03d.mat', P.zeroAngleIndex)));   % 0° 平面波、P1 体模
img = double(S.img);  w2 = double(S.w2);            % 波束形成信号与权重和

% ---- ROI：30 mm 附近 10 x 10 mm，即发射焦点附近 ------------------------
[XX, ZZ] = meshgrid(P.imageX, P.imageZ);
roi = abs(XX) <= 5e-3 & ZZ >= 25e-3 & ZZ <= 35e-3;

Ps = mean(img(roi).^2);                             % ROI 内平均信号功率
W  = mean(w2(roi));                                 % 平均权重平方和
% 要求 10*log10(Ps / (sigma_e^2 W)) = pwSnrTargetDB，反解 sigma_e
sigma_e = sqrt(Ps / (10^(P.pwSnrTargetDB/10) * W));
fprintf('signal power %.3e, mean sum(w^2) %.1f -> sigma_e %.3e\n', Ps, W, sigma_e);

% ---------------------------------------------------- 蒙特卡洛验证
% 用小网格是为了让 20 次重复足够便宜。每次生成独立的白色噪声场，送入同一个
% DAS，再把经验标准差与 sigma_e*sqrt(w2) 比较。
xr = (-5:0.1711:5)*1e-3;  zr = (25:0.1:35)*1e-3;
Nt = size(S.ch, 1);
% 传全零数据进去只是为了拿到 w2——权重与数据无关，这样省掉一次真实波束形成。
[~, w2r] = das_pw(zeros(Nt,P.nElements,'single'), S.t0, 0, P, xr, zr, S.dt);  % 只要权重
acc = zeros(numel(zr), numel(xr));  acc2 = acc;
rng(11,'twister');                                  % 固定种子，验证可复现
for r = 1:P.nMonteCarlo
    nz = sigma_e * randn(Nt, P.nElements);                  % 白色通道噪声
    nr = das_pw(nz, S.t0, 0, P, xr, zr, S.dt);      % 波束形成后的噪声
    acc  = acc  + nr;
    acc2 = acc2 + nr.^2;
end
% 经验标准差用 E[x^2]-E[x]^2 而不是 std()：只有 20 次实现，先减掉样本均值才
% 不会把均值的偏差算进方差；max(...,0) 防止数值误差把方差压成负数再开方。
emp = sqrt(max(acc2/P.nMonteCarlo - (acc/P.nMonteCarlo).^2, 0));    % 经验标准差（E[x^2]-E[x]^2）
ana = sigma_e * sqrt(w2r);                          % 解析标准差
ratio = emp ./ ana;
fprintf('Monte-Carlo (%d reps): empirical/analytic noise ratio = %.3f +- %.3f\n', ...
    P.nMonteCarlo, mean(ratio(:)), std(ratio(:)));

% ---- 追加写回 calib.mat（其余字段保持不变） ----------------------------
% Ps30 / W30 一起存：报告里若要复核 sigma_e，只要有这两个数和目标 SNR 就能重算。
C.sigma_e = sigma_e;
C.mcRatio = mean(ratio(:));
C.mcRatioStd = std(ratio(:));
C.Ps30 = Ps;  C.W30 = W;
save(fullfile(P.dataDir,'calib.mat'), '-struct', 'C');
fprintf('sigma_e saved to calib.mat\n');
end
