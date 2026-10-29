function stage_analysis(P)
%STAGE_ANALYSIS  E3/E4/E5：相干复合、指标计算与汇总表。
%
%   STAGE_ANALYSIS(P) 是整个实验的“规约（reduce）”步骤。它读取采集阶段写下的
%   逐次发射波束形成数据，输出：
%     data/res/results.mat  ：所有图像与所有指标（供画图使用）
%     data/res/*.csv        ：同样内容的 CSV（供生成报告表格）
%
%   输入（都在磁盘上）
%     data/pw/<TAG>_a###.mat   波束形成射频，每个角度/体模一个文件
%     data/foc/<TAG>_l###.mat  波束形成射频，每条聚焦线一个文件
%     data/noise/N1_a###.mat   波束形成后的噪声，每个角度一个文件
%
%   每个复合档位的处理链
%     1. 把该子集的 n 幅波束形成射频相加            （相干！）
%     2. 相加之后才取 Hilbert 包络                  （绝不能提前）
%     3. 对数压缩到 60 dB 用于显示
%     4. 测分辨率 / 对比度 / SNR
%
%   如果先取包络再相加，就是非相干复合——那正是论文反复强调的错误做法，横向
%   分辨率一点都不会改善。把求和留在射频域，是本函数最核心的实现细节。
%
%   指标约定（以及读表时要注意的地方）
%   ----------------------------------
%   * 分辨率只在点目标体模（P3）上测**横向** -6 dB 宽度。这个判据与 pulse_stats
%     的 -6 dB 带宽同口径（半幅度）；轴向 -6 dB 宽度虽然也算得出来，但三种方法
%     一样，不作为比较项（见下面“分辨率”一段）。
%   * 对比度用“强度”（包络平方），内区是半径 2 mm 的圆盘，外区是 4–6 mm 的环带。
%     无噪声图像上它只由旁瓣泄漏决定，实测 |C| = 12.8 / 23.3 / 16.0 dB（PW /
%     C11 / 单焦点），与论文的 12 / 20 / 12–32 dB 同量级；含噪版本（C_noisy）
%     另算一套，两者回答的问题不同（见下面“原理”一段）。
%   * CNR 用同样的区域；囊肿一黑，它就饱和在 0.7 附近，在本实验条件下并不能区分
%     方法优劣。列出它是为了完整，而不是当作排序依据。
%   * SNR = 10log10(<s^2>/<n^2>)，5x5 mm 滑窗，其中噪声功率来自真正经过同一个
%     DAS 的噪声实现（不是公式推的）。报**绝对** dB：刻度锚点由 stage_sigma 定，
%     σ_e 使单平面波在 30 mm 处为 20 dB（实测 20.4 dB），所以这一列数本身不是
%     临床可比的灵敏度；有意义的量是方法之间的差（C11−PW = 8.8 dB）。
%
%   数据依赖关系（哪一组输入决定哪一组输出）
%   ----------------------------------------
%     P1  的 11 个角度  ──►  PW 与 C11 的 B 模式、SNR
%     P3  的 11 个角度  ──►  横向分辨率（轴向宽度只作 QA，不作比较）
%     REM 的 11 个角度  ──►  与 P1 相减得到 P2，再算对比度/CNR
%     noise 的 11 帧     ──►  SNR 的噪声项、含噪对比度
%     foc/*.mat         ──►  单焦点参照（分辨率、对比度、SNR 各一条曲线）
%   注意 P2 不是独立采集出来的，而是 P1 − REM；这一减法只在射频域、取包络之前
%   做才严格成立。
%
%   原理：复合为什么能提高 SNR
%   -------------------------
%   DAS 已经把每个像素的几何延时补掉，所以同一个散射点的回波在 11 个角度下是
%   同相的：相加时信号幅度变 11 倍（功率 121 倍）。而每次发射的噪声彼此独立：
%   相加时功率只变 11 倍。于是 SNR 从 S/N 变成 11S/(sqrt(11)N) = sqrt(11)·S/N，
%   即 10log10(11) = 10.4 dB。本实验实测 8.8 dB，差的那 1.6 dB 来自各角度的信号
%   幅度并不完全相等（偏转损失、孔径权重差异）与噪声并非严格白噪声。这张"理论上
%   限"就写在 R.snrGainTheory 里，报告用它作对照。
%
%   原理：为什么 P2 = P1 − REM 必须在射频域做
%   -----------------------------------------
%   取包络是取模，是非线性的：|a| − |b| 与 |a − b| 不是一回事。P1 与 REM 两幅图
%   的散斑相位不同，模值相减得到的东西没有物理意义；只有在射频域相减，才等价于
%   "把囊肿内那批散射子的贡献从回波里去掉"。
%
%   原理：为什么对比度要算两套（无噪 / 含噪）
%   ----------------------------------------
%   无噪对比度衡量旁瓣泄漏：囊肿内部本该是零，实际有值完全是旁瓣把外面的强回波
%   漏了进来，所以它反映的是波束形成质量。含噪对比度衡量噪声底：真实系统里囊肿
%   内部的"黑"由噪声决定，所以它反映的是成像灵敏度。两者都可能成为限制因素，
%   只报一个会得出片面结论。
%
%   原理：CNR 为什么对本实验不敏感
%   ------------------------------
%   CNR = (μ_out − μ_in)/σ_out。囊肿一黑（μ_in ≈ 0），它就趋近 μ_out/σ_out ——
%   对完全发展的散斑，这个比值近似是常数（本实验 0.54–0.63），与用了什么方法关系
%   很小。所以 CNR 列出来只是为了完整，不作为方法排序的依据。

fprintf('\n===== analysis =====\n');
outdir = fullfile(P.dataDir, 'res');
if ~exist(outdir,'dir'), mkdir(outdir); end

idxAll = subset_indices(P);          % 两种方法各自的角度下标（1和[1, 2, ..., 11]）
% 名字按 P.compoundCounts 自动生成：n=1 叫 PW，n>1 叫 Cn；顺序与 P.compoundCounts 一致。
names = cell(1, numel(P.compoundCounts));
for k = 1:numel(P.compoundCounts)
    if P.compoundCounts(k) == 1
        names{k} = 'PW';
    else
        names{k} = sprintf('C%d', P.compoundCounts(k));
    end
end
nM = numel(names);

R = struct();
R.names = names;
R.nAngles = P.compoundCounts;
R.x = P.imageX;  R.z = P.imageZ;

% ---------------------------------------------------------------- 复合 P1
% 每档一遍：用同一组角度下标同时复合三个体模与噪声，这样 P2 = P1 - REM 自洽，
% 噪声也正好属于这几次发射。
for k = 1:nM
    rf1 = compound_rf(P, 'P1',  idxAll{k});      % 均匀散斑体模
    rfR = compound_rf(P, 'REM', idxAll{k});      % 囊肿内部被移除的散射子
    rf3 = compound_rf(P, 'P3',  idxAll{k});      % 点目标
    nz1 = compound_noise(P, idxAll{k});          % 同样这几次发射的噪声

    % ---- 包络必须在复合之后取 ------------------------------------------
    env1 = abs(hilbert(double(rf1)));
    env3 = abs(hilbert(double(rf3)));
    env2 = abs(hilbert(double(rf1 - rfR)));      % P2 = P1 - REM（幅度线性）
    env2n= abs(hilbert(double(rf1 - rfR + nz1)));% 同上，叠加采集噪声

    R.P1.rf(:,:,k) = rf1;                        % 原始射频，图 F7 需要
    R.P1.nz(:,:,k) = nz1;                        % 噪声射频，SNR 指标需要
    R.P1.db(:,:,k) = to_db(env1, P.displayRangeDB);        % B 模式，60 dB
    R.P2.db(:,:,k) = to_db(env2, P.displayRangeDB);
    R.P3.db(:,:,k) = to_db(env3, P.displayRangeDB);
    R.P3.env(:,:,k) = single(env3);              % 包络留着算 PSF 指标

    R.contrast(k)      = metrics_contrast(env2,  P.imageX, P.imageZ, P);   % 无噪声
    R.contrastN(k)     = metrics_contrast(env2n, P.imageX, P.imageZ, P);   % 含噪声

    fprintf('  %s done\n', names{k});
end

% ------------------------------------------------------------ 聚焦图像
% 传统 B 模式按线拼接成二维图像。P3 用 256 线（PSF 采样更细），P1/P2 用 128 线
% （一个阵元一条线）。
rfF1 = load_foc(P, 'P1',  P.txLineX);
rfFR = load_foc(P, 'REM', P.txLineX);
rfF3 = load_foc(P, 'P3',  P.psfLineX);
nzF  = load_foc_noise(P, P.txLineX);

R.FOC.x = P.txLineX;
R.FOC.db1   = to_db(abs(hilbert(double(rfF1))), P.displayRangeDB);
R.FOC.db2   = to_db(abs(hilbert(double(rfF1 - rfFR))), P.displayRangeDB);
R.FOC.db2n  = to_db(abs(hilbert(double(rfF1 - rfFR + nzF))), P.displayRangeDB);
R.FOC.env3  = single(abs(hilbert(double(rfF3))));
% env3 是 256 线（P.psfLineX，0.171 mm）的重建结果，而上面的 B 模式图是 128 线
% （P.txLineX，0.33 mm）：两者的横向轴**不一样长**，所以必须各存各的。
% 这个字段是给 env3 配的轴；把 env3 和 R.FOC.x 一起用会把线间距错算一倍。
R.FOC.xP3   = P.psfLineX;
R.FOC.nz    = nzF;
R.FOC.rf1   = rfF1;

envF2  = abs(hilbert(double(rfF1 - rfFR)));
envF2n = abs(hilbert(double(rfF1 - rfFR + nzF)));
R.contrastF  = metrics_contrast(envF2,  P.txLineX, P.imageZ, P);
R.contrastFN = metrics_contrast(envF2n, P.txLineX, P.imageZ, P);

% ------------------------------------------------------------- 分辨率
% 轴线上六个目标（15…55 mm），得到“横向分辨率随深度”曲线。
%
% 只比横向。轴向分辨率由脉冲的空间长度 c*tau/2 决定，与聚焦方式和复合都无关：
% 三种方法实测都是 0.19 mm（0.19—0.20 mm，±3 %，落在读数不确定度以内），比不出
% 差别。metrics_resolution 仍然返回 wAx（QA 行用它检查 PSF 是不是约 0.2 mm），但
% 不作为比较项进入结果与表格——原因记在 docs/实验报告_讨论与附录.md 的 B14。
zPts = [15 25 30 35 45 55]*1e-3;
xPts = [0 0 0 0 0 0]*1e-3;
R.zPts = zPts;  R.xPts = xPts;
wLat = nan(nM, numel(zPts));
for k = 1:nM
    for q = 1:numel(zPts)
        m = metrics_resolution(R.P3.env(:,:,k), P.imageX, P.imageZ, xPts(q), zPts(q));
        wLat(k,q) = m.wLat;
    end
end
wLatF = nan(1, numel(zPts));
for q = 1:numel(zPts)
    m = metrics_resolution(R.FOC.env3, P.psfLineX, P.imageZ, xPts(q), zPts(q));
    wLatF(q) = m.wLat;
end
R.wLat = wLat;  R.wLatF = wLatF;

% z = 30 mm 处的三个横向位置：检查分辨率在横向上是否均匀
xl = [-12 0 12]*1e-3;
R.wLatX = nan(nM, 3);  R.wLatXF = nan(1,3);
for q = 1:3
    for k = 1:nM
        m = metrics_resolution(R.P3.env(:,:,k), P.imageX, P.imageZ, xl(q), 30e-3);
        R.wLatX(k,q) = m.wLat;
    end
    m = metrics_resolution(R.FOC.env3, P.psfLineX, P.imageZ, xl(q), 30e-3);
    R.wLatXF(q) = m.wLat;
end

% ------------------------------------------------------------------- SNR
% SNR(x,z) = 10log10( <s^2> / <n^2> )，5x5 mm 滑窗；沿 x 取中位数得到每个深度
% 一条曲线。三种方法都报**绝对** SNR（dB），不再减掉某个参考点。
%
% 绝对刻度的锚点：σ_e 由 stage_sigma 按"单平面波在 30 mm 处 20 dB"标定，所以
% 单平面波在 30 mm 处读到约 20 dB（实测 20.4，差 0.45 来自 5x5 mm 平滑与散斑
% 统计，见报告 5.7）。因此绝对值的物理含义是"以这个标定约定为刻度"，不是临床
% 可比的绝对灵敏度；方法之间的差（如 C11 比 PW 高 8.8 dB）才是有意义的量。
% snrAbsPW30 只作为刻度锚点保留，供图和正文引用，不再用来做归一化。
%
% 为什么可以直接用射频的平方：5 mm 窗内有十几个波长（λ = 0.37 mm），求平均时
% 载波振荡被平滑掉，<rf^2> 正比于包络功率，不必每帧再取一次 Hilbert。
% 为什么沿 x 取中位数而不是均值：散斑让单条线的 SNR 起伏很大，中位数不会被少数
% 强散射点拉偏；这也与"沿横线取一个代表值"的物理直觉一致。
zdep = (15:5:55)*1e-3;
[zz, snrPW] = snr_curve(R.P1.rf(:,:,1), R.P1.nz(:,:,1), P, P.pixelPitchX);
snrM = zeros(nM, numel(zdep));
for k = 1:nM
    [~, s] = snr_curve(R.P1.rf(:,:,k), R.P1.nz(:,:,k), P, P.pixelPitchX);
    snrM(k,:) = interp1(zz, s, zdep, 'linear');
end
[zzF, sF] = snr_curve(rfF1, nzF, P, P.txLinePitch);
snrF = interp1(zzF, sF, zdep, 'linear');
R.zdep = zdep;
R.snr  = snrM;                                        % 绝对 SNR，[nM x 9]，dB
R.snrF = snrF;                                        % 单焦点绝对 SNR，[1 x 9]
R.snrAbsPW30 = interp1(zz, snrPW, 30e-3, 'linear');   % 刻度锚点（标定目标 20 dB）

% ------------------------------------------------- 理论增益（用于前向检查）
% 理想相干复合增益；报告把实测值与之对比。
R.snrGainTheory = 10*log10(P.compoundCounts);

% ------------------------------------------------------------------ 保存
save(fullfile(outdir,'results.mat'), '-struct', 'R', '-v7.3');
fprintf('saved %s\n', fullfile(outdir,'results.mat'));

write_csv(outdir, names, R);         % 同样内容导出 CSV，供报告使用
end

% =========================================================================
function rf = compound_rf(P, tag, idx)
%COMPOUND_RF  若干角度的波束形成射频相干求和（单个体模）。
%   求和刻意放在射频域：对包络求和属于非相干复合，不会锐化发射波束。
rf = zeros(numel(P.imageZ), numel(P.imageX), 'single');
for i = idx
    S = load(fullfile(P.dataDir,'pw',sprintf('%s_a%03d.mat', tag, i)));
    rf = rf + S.img;
end
end

function rf = compound_noise(P, idx)
%COMPOUND_NOISE  同样这些角度的噪声图像相干求和。
rf = zeros(numel(P.imageZ), numel(P.imageX), 'single');
for i = idx
    S = load(fullfile(P.dataDir,'noise',sprintf('N1_a%03d.mat', i)));
    rf = rf + S.img;
end
end

function rf = load_foc(P, tag, xf)
%LOAD_FOC  把某个体模的逐线文件拼成完整聚焦 B 模式图像。
rf = zeros(numel(P.imageZ), numel(xf), 'single');
for k = 1:numel(xf)
    S = load(fullfile(P.dataDir,'foc',sprintf('%s_l%03d.mat', tag, k)));
    rf(:,k) = S.img;
end
end

function rf = load_foc_noise(P, xf)
%LOAD_FOC_NOISE  同上，拼接只有噪声的聚焦图像。
rf = zeros(numel(P.imageZ), numel(xf), 'single');
for k = 1:numel(xf)
    S = load(fullfile(P.dataDir,'noise',sprintf('NF_l%03d.mat', k)));
    rf(:,k) = S.img;
end
end

function [zz, snr] = snr_curve(rf, nrf, P, dx)
%SNR_CURVE  逐像素 SNR 的深度曲线（dB）。
%   信号与噪声功率都用 5x5 mm 盒平滑（可分离，所以用两个一维滤波器而不是一个
%   二维卷积）。沿 x 取中位数：散斑让单个像素的 SNR 起伏很大，中位数更稳。
ps = double(rf).^2;  pn = double(nrf).^2;

h1 = ones(round(5e-3/P.pixelPitchZ),1); h1 = h1/sum(h1);     % 轴向 5 mm 盒
h2 = ones(1, round(5e-3/dx));  h2 = h2/sum(h2);     % 横向 5 mm 盒
ps = imfilter(imfilter(ps, h1, 'replicate'), h2, 'replicate');
pn = imfilter(imfilter(pn, h1, 'replicate'), h2, 'replicate');

S  = 10*log10(ps ./ max(pn, realmin));               % SNR 图 [dB]
zz = P.imageZ;
snr = median(S, 2);                                  % 每个深度一个值
end

function write_csv(outdir, names, R)
%WRITE_CSV  把指标表导出为 CSV，这样重新生成报告时不必再跑 MATLAB 分析。
%   每类指标一个文件。

% ---- 横向分辨率：六个轴上深度处的 -6 dB 宽度 --------------------------
fid = fopen(fullfile(outdir,'table_resolution.csv'),'w');
fprintf(fid, 'method,');
for q = 1:numel(R.zPts), fprintf(fid,'wLat_%dmm,', round(R.zPts(q)*1e3)); end
fprintf(fid,'\n');
for k = 1:numel(names)
    fprintf(fid,'%s,', names{k});
    fprintf(fid,'%.3f,', R.wLat(k,:)*1e3);      % 单位 mm
    fprintf(fid,'\n');
end
fprintf(fid,'FOCUSED,');
fprintf(fid,'%.3f,', R.wLatF*1e3);
fprintf(fid,'\n');
fclose(fid);

% ---- 对比度与 CNR，无噪声与含噪声各一份 ------------------------------
fid = fopen(fullfile(outdir,'table_contrast.csv'),'w');
fprintf(fid,'method,nAngles,C_dB,C_std,CNR,CNR_std,C_noisy_dB,CNR_noisy\n');
for k = 1:numel(names)
    fprintf(fid,'%s,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', names{k}, R.nAngles(k), ...
        R.contrast(k).Cmean, R.contrast(k).Cstd, R.contrast(k).CNRmean, R.contrast(k).CNRstd, ...
        R.contrastN(k).Cmean, R.contrastN(k).CNRmean);
end
fprintf(fid,'FOCUSED,128,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', ...
    R.contrastF.Cmean, R.contrastF.Cstd, R.contrastF.CNRmean, R.contrastF.CNRstd, ...
    R.contrastFN.Cmean, R.contrastFN.CNRmean);
fclose(fid);

% ---- 绝对 SNR 随深度 ------------------------------------------------
fid = fopen(fullfile(outdir,'table_snr.csv'),'w');
fprintf(fid,'method,'); fprintf(fid,'%dmm,', round(R.zdep*1e3)); fprintf(fid,'\n');
for k = 1:numel(names)
    fprintf(fid,'%s,', names{k});
    fprintf(fid,'%.2f,', R.snr(k,:));
    fprintf(fid,'\n');
end
fprintf(fid,'FOCUSED,'); fprintf(fid,'%.2f,', R.snrF); fprintf(fid,'\n');
fclose(fid);
end
