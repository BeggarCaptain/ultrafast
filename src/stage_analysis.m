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
%   * 分辨率只在点目标体模（P3）上测**横向** -10 dB 宽度。这个电平是设计文档
%     确定的；用混了（例如取 -6 dB）会让数值差 30 % 左右。轴向 -6 dB 宽度虽然
%     也算得出来，但三种方法一样，不作为比较项（见下面“分辨率”一段）。
%   * 对比度用“强度”（包络平方），内区是半径 2 mm 的圆盘，外区是 4–6 mm 的环带。
%     在无噪声图像上它由旁瓣决定，所以这里能到 40 dB 以上（11 角度实测 45 dB），
%     而论文只有 12–33 dB（论文的对比度是噪声受限的）。
%   * CNR 用同样的区域；囊肿一黑，它就饱和在 0.7 附近，在本实验条件下并不能区分
%     方法优劣。列出它是为了完整，而不是当作排序依据。
%   * SNR = 10log10(<s^2>/<n^2>)，5x5 mm 滑窗，其中噪声功率来自真正经过同一个
%     DAS 的噪声实现（不是公式推的）。全部相对于 30 mm 处的单平面波，这也是论文
%     报告复合增益的方式。
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
% 差别。metrics_resolution 仍然返回 w6（QA 行用它检查 PSF 是不是约 0.2 mm），但
% 不作为比较项进入结果与表格——原因记在 docs/实验报告_讨论与附录.md 的 B14。
zPts = [15 25 30 35 45 55]*1e-3;
xPts = [0 0 0 0 0 0]*1e-3;
R.zPts = zPts;  R.xPts = xPts;
w10 = nan(nM, numel(zPts));
for k = 1:nM
    for q = 1:numel(zPts)
        m = metrics_resolution(R.P3.env(:,:,k), P.imageX, P.imageZ, xPts(q), zPts(q));
        w10(k,q) = m.w10;
    end
end
w10F = nan(1, numel(zPts));
for q = 1:numel(zPts)
    m = metrics_resolution(R.FOC.env3, P.psfLineX, P.imageZ, xPts(q), zPts(q));
    w10F(q) = m.w10;
end
R.w10 = w10;  R.w10F = w10F;

% z = 30 mm 处的三个横向位置：检查分辨率在横向上是否均匀
xl = [-12 0 12]*1e-3;
R.w10Lat = nan(nM, 3);  R.w10LatF = nan(1,3);
for q = 1:3
    for k = 1:nM
        m = metrics_resolution(R.P3.env(:,:,k), P.imageX, P.imageZ, xl(q), 30e-3);
        R.w10Lat(k,q) = m.w10;
    end
    m = metrics_resolution(R.FOC.env3, P.psfLineX, P.imageZ, xl(q), 30e-3);
    R.w10LatF(q) = m.w10;
end

% ------------------------------------------------------------------- SNR
% SNR(x,z) = 10log10( <s^2> / <n^2> )，5x5 mm 滑窗；沿 x 取中位数得到每个深度
% 一条曲线。最后全部相对于 30 mm 处的单平面波——这就是论文报告增益的方式。
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
R.snr  = snrM;
R.snrF = snrF;
R.snrAbsPW30 = interp1(zz, snrPW, 30e-3, 'linear');   % 0 dB 参考点
R.snrRel = snrM - R.snrAbsPW30;
R.snrRelF = snrF - R.snrAbsPW30;

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

% ---- 横向分辨率：六个轴上深度处的 -10 dB 宽度 -------------------------
fid = fopen(fullfile(outdir,'table_resolution.csv'),'w');
fprintf(fid, 'method,');
for q = 1:numel(R.zPts), fprintf(fid,'w10_%dmm,', round(R.zPts(q)*1e3)); end
fprintf(fid,'\n');
for k = 1:numel(names)
    fprintf(fid,'%s,', names{k});
    fprintf(fid,'%.3f,', R.w10(k,:)*1e3);      % 单位 mm
    fprintf(fid,'\n');
end
fprintf(fid,'FOCUSED,');
fprintf(fid,'%.3f,', R.w10F*1e3);
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

% ---- 相对 SNR 随深度 ------------------------------------------------
fid = fopen(fullfile(outdir,'table_snr.csv'),'w');
fprintf(fid,'method,'); fprintf(fid,'%dmm,', round(R.zdep*1e3)); fprintf(fid,'\n');
for k = 1:numel(names)
    fprintf(fid,'%s,', names{k});
    fprintf(fid,'%.2f,', R.snrRel(k,:));
    fprintf(fid,'\n');
end
fprintf(fid,'FOCUSED,'); fprintf(fid,'%.2f,', R.snrRelF); fprintf(fid,'\n');
fclose(fid);
end
