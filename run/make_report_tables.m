function make_report_tables()
%MAKE_REPORT_TABLES  把报告要用的表格导出为 Markdown（data/res/tables.md）。
%
%   报告里绝不允许出现手敲的数字：每一个数值都来自这个文件，而它直接由
%   results.mat 与 calib.mat 生成。指标一变，重跑本函数，表格就重新一致了。
%
%   生成的表
%     T-A ：E0 标定与自检
%     T-B ：横向分辨率（-10 dB）随深度
%     T-C ：轴向分辨率（-6 dB）随深度
%     T-D ：囊肿对比度与 CNR，无噪声与含噪声
%     T-E ：相对 SNR 随深度，附 10log10(n) 参考
%     T-F ：散斑统计（强度与包络的 μ/σ、相关长度）
%
%   文件以 UTF-8 打开，中文表头与 ±、µ 等符号才不会乱码；onCleanup 保证即使
%   中途出错也会关闭文件句柄。
%
%   表格口径（与报告正文一致，读表时请对照）
%   ----------------------------------------
%   * T-B / T-C 的列是点目标所在的六个深度；最后一行“单焦点”是常规逐线聚焦，
%     不是复合档位，所以它没有“n 角度”的概念。
%   * T-D 的 |C| 取绝对值（原始 C 为负）；“含噪”指叠加了采集噪声的图像，两者
%     共用同一份噪声实现，可比。
%   * T-E 全部以“单平面波在 30 mm 处”为 0 dB；最后一行 10log₁₀(n) 与深度无关，
%     所以单独写成一行而不是塞进表里。
%   * T-F 同时给强度与包络两个 μ/σ：前者理想值 1.00（指数分布），后者 1.91
%     （Rayleigh）。只报一个数最容易引起误会。

here = fileparts(mfilename('fullpath'));   % .../ultrafast/run
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);
C = load(fullfile(P.dataDir, 'calib.mat'));
R = load(fullfile(P.dataDir, 'res', 'results.mat'));

fid = fopen(fullfile(P.dataDir, 'res', 'tables.md'), 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));      % 即使下面某行 fprintf 出错也会关闭

nM = numel(R.names);

% ---------------------------------------------------------------- T-A：E0
fprintf(fid, '### T-A 标定与自检（E0）\n\n');
fprintf(fid, '| 项目 | 结果 | 判据 |\n|---|---|---|\n');
fprintf(fid, '| 中心频率 | %.2f MHz | 4.5 ± 0.1 |\n', C.pulse.fc/1e6);
fprintf(fid, '| −6 dB 相对带宽 | %.1f %% | 60 ± 5 |\n', 100*C.pulse.bw);
fprintf(fid, '| 倾角斜率比（%g°） | %.3f | 1 ± 0.05 |\n', C.alphaChk(2)*180/pi, C.tiltRatio(2));
fprintf(fid, '| 倾角斜率比（%g°） | %.3f | 1 ± 0.05 |\n', C.alphaChk(3)*180/pi, C.tiltRatio(3));
fprintf(fid, '| 时间零点校正 dt | %.3f–%.3f µs（极差 %.3f） | 常数 |\n', ...
    min(C.dt)*1e6, max(C.dt)*1e6, (max(C.dt)-min(C.dt))*1e6);
fprintf(fid, '| 线性/样条插值差 | %.1f dB | < 0.5 dB 量级 |\n', C.interpErrDB);
fprintf(fid, '| 20k 散射子单次发射 | %.1f s | — |\n', C.tAcq20k);
fprintf(fid, '| 噪声标定 σ_e | %.3e | PW@30 mm = %d dB |\n', C.sigma_e, P.snrPW30);
fprintf(fid, '| 噪声解析/蒙特卡洛比 | %.3f ± %.3f | 1.0 |\n', C.mcRatio, C.mcRatioStd);

% ------------------------------------------------- T-B：横向分辨率
% 每个方法一行，每个点目标深度一列；最后一行是常规聚焦图像（256 线），单独列出。
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

% --------------------------------------------------- T-C：轴向分辨率
% 同样的布局，-6 dB 宽度。它应当与方法无关：复合改变的是横向发射孔径，
% 并不会改变脉冲的轴向带宽。
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

% ------------------------------------------------------ T-D：对比度 / CNR
% 九个囊肿的均值 ± 标准差，无噪声与含噪声各一组。
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

% ------------------------------------------------------------- T-E：SNR
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
fprintf(fid, '\n');
% 理论值 10*log10(n) 与深度无关，所以单独写成一小行，而不是塞进按深度排的表里
fprintf(fid, '\n10log10(n) 理论：');
for k = 1:nM
    fprintf(fid, ' %s=%.1f dB；', R.names{k}, 10*log10(R.nAngles(k)));
end
fprintf(fid, '\n');

% --------------------------------------------------------- T-F：散斑
% 两个比值是有意并排的：1.00 是“强度”的理想值（指数分布），1.91 是“包络”的
% 理想值（Rayleigh）。
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
