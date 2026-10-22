function make_report_tables()
%MAKE_REPORT_TABLES  Emit markdown tables (data/res/tables.md) from the
%   numeric results, so that the report text never has to re-type numbers.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);
C = load(fullfile(P.dataDir, 'calib.mat'));
R = load(fullfile(P.dataDir, 'res', 'results.mat'));

fid = fopen(fullfile(P.dataDir, 'res', 'tables.md'), 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));

names = [R.names, {'单焦点'}];
nM = numel(R.names);

fprintf(fid, '### T-A 标定与自检（E0）\n\n');
fprintf(fid, '| 项目 | 结果 | 判据 |\n|---|---|---|\n');
fprintf(fid, '| 中心频率 | %.2f MHz | 4.5 ± 0.1 |\n', C.pulse.fc/1e6);
fprintf(fid, '| −6 dB 相对带宽 | %.1f %% | 60 ± 5 |\n', 100*C.pulse.bw);
fprintf(fid, '| 倾角斜率比（10°） | %.3f | 1 ± 0.05 |\n', C.tiltRatio(2));
fprintf(fid, '| 倾角斜率比（16.5°） | %.3f | 1 ± 0.05 |\n', C.tiltRatio(3));
fprintf(fid, '| 时间零点校正 dt | %.3f–%.3f µs（极差 %.3f） | 常数 |\n', ...
    min(C.dt)*1e6, max(C.dt)*1e6, (max(C.dt)-min(C.dt))*1e6);
fprintf(fid, '| 线性/样条插值差 | %.1f dB | < 0.5 dB 量级 |\n', C.interpErrDB);
fprintf(fid, '| 20k 散射子单次发射 | %.1f s | — |\n', C.tAcq20k);
fprintf(fid, '| 噪声标定 σ_e | %.3e | PW@30 mm = %d dB |\n', C.sigma_e, P.snrPW30);
fprintf(fid, '| 噪声解析/蒙特卡洛比 | %.3f ± %.3f | 1.0 |\n', C.mcRatio, C.mcRatioStd);

fprintf(fid, '\n### T-B 横向分辨率（−10 dB 宽度，mm）\n\n');
fprintf(fid, '| 方法 |');
fprintf(fid, ' %d mm |', round(R.zPts*1e3));
fprintf(fid, '\n|---|'); fprintf(fid, '---|'); fprintf(fid, '\n');
for k = 1:nM
    fprintf(fid, '| %s |', R.names{k});
    fprintf(fid, ' %.2f |', R.w10(k,:)*1e3);
    fprintf(fid, '\n');
end
fprintf(fid, '| 单焦点 |');
fprintf(fid, ' %.2f |', R.w10F*1e3);
fprintf(fid, '\n');

fprintf(fid, '\n### T-C 轴向分辨率（−6 dB 宽度，mm）\n\n');
fprintf(fid, '| 方法 |');
fprintf(fid, ' %d mm |', round(R.zPts*1e3));
fprintf(fid, '\n|---|'); fprintf(fid, '---|'); fprintf(fid, '\n');
for k = 1:nM
    fprintf(fid, '| %s |', R.names{k});
    fprintf(fid, ' %.2f |', R.w6(k,:)*1e3);
    fprintf(fid, '\n');
end
fprintf(fid, '| 单焦点 |');
fprintf(fid, ' %.2f |', R.w6F*1e3);
fprintf(fid, '\n');

fprintf(fid, '\n### T-D 对比度与 CNR（9 个囊肿，均值±标准差）\n\n');
fprintf(fid, '| 方法 | n | |C| 无噪 (dB) | |C| 含噪 (dB) | CNR 无噪 | CNR 含噪 |\n|---|---|---|---|---|---|\n');
for k = 1:nM
    fprintf(fid, '| %s | %d | %.1f ± %.1f | %.1f ± %.1f | %.2f ± %.2f | %.2f ± %.2f |\n', ...
        R.names{k}, R.nAngles(k), R.contrast(k).Cmean, R.contrast(k).Cstd, ...
        R.contrastN(k).Cmean, R.contrastN(k).Cstd, ...
        R.contrast(k).CNRmean, R.contrast(k).CNRstd, ...
        R.contrastN(k).CNRmean, R.contrastN(k).CNRstd);
end
fprintf(fid, '| 单焦点 | 128 | %.1f ± %.1f | %.1f ± %.1f | %.2f ± %.2f | %.2f ± %.2f |\n', ...
    R.contrastF.Cmean, R.contrastF.Cstd, R.contrastFN.Cmean, R.contrastFN.Cstd, ...
    R.contrastF.CNRmean, R.contrastF.CNRstd, R.contrastFN.CNRmean, R.contrastFN.CNRstd);

fprintf(fid, '\n### T-E 相对 SNR（dB，以 PW 在 30 mm 为 0 dB）\n\n');
fprintf(fid, '| 方法 |');
fprintf(fid, ' %d mm |', round(R.zdep*1e3));
fprintf(fid, '\n|---|'); fprintf(fid, '---|'); fprintf(fid, '\n');
for k = 1:nM
    fprintf(fid, '| %s |', R.names{k});
    fprintf(fid, ' %.1f |', R.snrRel(k,:));
    fprintf(fid, '\n');
end
fprintf(fid, '| 单焦点 |');
fprintf(fid, ' %.1f |', R.snrRelF);
fprintf(fid, '\n| 10log10(n) 理论 |');
fprintf(fid, ' %.1f |', 10*log10(R.nAngles));
fprintf(fid, '\n');

fprintf(fid, '\n### T-F 散斑统计（10×10 mm ROI，z = 30 mm）\n\n');
fprintf(fid, '| 方法 | 强度 μ/σ | 包络 μ/σ | 横向相关长度 (mm) | 轴向相关长度 (mm) |\n|---|---|---|---|---|\n');
for k = 1:nM
    fprintf(fid, '| %s | %.2f | %.2f | %.2f | %.2f |\n', R.names{k}, R.speck(k).ratio, ...
        R.speck(k).ratioA, R.speck(k).corrX*1e3, R.speck(k).corrZ*1e3);
end
fprintf(fid, '| 单焦点 | %.2f | %.2f | %.2f | %.2f |\n', R.speckF.ratio, R.speckF.ratioA, ...
    R.speckF.corrX*1e3, R.speckF.corrZ*1e3);
fprintf(fid, '\n（完全发展散斑：强度 μ/σ = 1.00（指数分布），包络 μ/σ = 1.91（Rayleigh））\n');
end
