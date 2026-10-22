function stage_analysis(P)
%STAGE_ANALYSIS  E3/E4/E5: compounding, metrics, tables.
%   Reads data/pw, data/foc, data/noise; writes data/res/results.mat + CSV.

fprintf('\n===== analysis =====\n');
outdir = fullfile(P.dataDir, 'res');
if ~exist(outdir,'dir'), mkdir(outdir); end

idxAll = subset_indices(P);
names = {'PW','C5','C9','C17','C33','C45','C71'};
nM = numel(names);

R = struct();
R.names = names;
R.nAngles = P.subsetN;
R.x = P.x;  R.z = P.z;

% ---------------------------------------------------------------- compound P1
for k = 1:nM
    rf1 = compound_rf(P, 'P1',  idxAll{k});
    rfR = compound_rf(P, 'REM', idxAll{k});
    rf3 = compound_rf(P, 'P3',  idxAll{k});
    nz1 = compound_noise(P, idxAll{k});

    env1 = abs(hilbert(double(rf1)));
    env3 = abs(hilbert(double(rf3)));
    env2 = abs(hilbert(double(rf1 - rfR)));
    env2n= abs(hilbert(double(rf1 - rfR + nz1)));

    R.P1.rf(:,:,k) = rf1;
    R.P1.nz(:,:,k) = nz1;
    R.P1.db(:,:,k) = to_db(env1, P.drDB);
    R.P2.db(:,:,k) = to_db(env2, P.drDB);
    R.P3.db(:,:,k) = to_db(env3, P.drDB);
    R.P3.env(:,:,k) = single(env3);

    R.contrast(k)      = metrics_contrast(env2,  P.x, P.z, P);
    R.contrastN(k)     = metrics_contrast(env2n, P.x, P.z, P);
    R.speck(k)         = metrics_speckle(env1, P.x, P.z, 0, 30e-3, 5e-3);

    fprintf('  %s done\n', names{k});
end

% ------------------------------------------------------------ focused images
rfF1 = load_foc(P, 'P1',  P.xf);
rfFR = load_foc(P, 'REM', P.xf);
rfF3 = load_foc(P, 'P3',  P.xpsf);
nzF  = load_foc_noise(P, P.xf);

R.FOC.x = P.xf;
R.FOC.db1   = to_db(abs(hilbert(double(rfF1))), P.drDB);
R.FOC.db2   = to_db(abs(hilbert(double(rfF1 - rfFR))), P.drDB);
R.FOC.db2n  = to_db(abs(hilbert(double(rfF1 - rfFR + nzF))), P.drDB);
R.FOC.env3  = single(abs(hilbert(double(rfF3))));
R.FOC.nz    = nzF;
R.FOC.rf1   = rfF1;

envF2  = abs(hilbert(double(rfF1 - rfFR)));
envF2n = abs(hilbert(double(rfF1 - rfFR + nzF)));
R.contrastF  = metrics_contrast(envF2,  P.xf, P.z, P);
R.contrastFN = metrics_contrast(envF2n, P.xf, P.z, P);
R.speckF     = metrics_speckle(abs(hilbert(double(rfF1))), P.xf, P.z, 0, 30e-3, 5e-3);

% ------------------------------------------------------------- resolution
zPts = [15 25 30 35 45 55]*1e-3;
xPts = [0 0 0 0 0 0]*1e-3;
R.zPts = zPts;  R.xPts = xPts;
w10 = nan(nM, numel(zPts));  w6 = nan(nM, numel(zPts));
for k = 1:nM
    for q = 1:numel(zPts)
        m = metrics_resolution(R.P3.env(:,:,k), P.x, P.z, xPts(q), zPts(q));
        w10(k,q) = m.w10;  w6(k,q) = m.w6;
    end
end
w10F = nan(1, numel(zPts));  w6F = nan(1, numel(zPts));
for q = 1:numel(zPts)
    m = metrics_resolution(R.FOC.env3, P.xpsf, P.z, xPts(q), zPts(q));
    w10F(q) = m.w10;  w6F(q) = m.w6;
end
R.w10 = w10;  R.w6 = w6;  R.w10F = w10F;  R.w6F = w6F;

% lateral points at z = 30 mm
xl = [-12 0 12]*1e-3;
R.w10Lat = nan(nM, 3);  R.w10LatF = nan(1,3);
for q = 1:3
    for k = 1:nM
        m = metrics_resolution(R.P3.env(:,:,k), P.x, P.z, xl(q), 30e-3);
        R.w10Lat(k,q) = m.w10;
    end
    m = metrics_resolution(R.FOC.env3, P.xpsf, P.z, xl(q), 30e-3);
    R.w10LatF(q) = m.w10;
end

% ------------------------------------------------------------------- SNR
zdep = (15:5:55)*1e-3;
[zz, snrPW] = snr_curve(R.P1.rf(:,:,1), R.P1.nz(:,:,1), P, P.dx);
snrM = zeros(nM, numel(zdep));
for k = 1:nM
    [~, s] = snr_curve(R.P1.rf(:,:,k), R.P1.nz(:,:,k), P, P.dx);
    snrM(k,:) = interp1(zz, s, zdep, 'linear');
end
[zzF, sF] = snr_curve(rfF1, nzF, P, P.txPitch);
snrF = interp1(zzF, sF, zdep, 'linear');
R.zdep = zdep;
R.snr  = snrM;
R.snrF = snrF;
R.snrAbsPW30 = interp1(zz, snrPW, 30e-3, 'linear');
R.snrRel = snrM - R.snrAbsPW30;
R.snrRelF = snrF - R.snrAbsPW30;

% ------------------------------------------------- forward model check (E0/7)
R.snrGainTheory = 10*log10(P.subsetN);

% ------------------------------------------------------------------ save
R.hist = 'placeholder';
save(fullfile(outdir,'results.mat'), '-struct', 'R', '-v7.3');
fprintf('saved %s\n', fullfile(outdir,'results.mat'));

write_csv(outdir, names, R);
end

% =========================================================================
function rf = compound_rf(P, tag, idx)
rf = zeros(numel(P.z), numel(P.x), 'single');
for i = idx
    S = load(fullfile(P.dataDir,'pw',sprintf('%s_a%03d.mat', tag, i)));
    rf = rf + S.img;
end
end

function rf = compound_noise(P, idx)
rf = zeros(numel(P.z), numel(P.x), 'single');
for i = idx
    S = load(fullfile(P.dataDir,'noise',sprintf('N1_a%03d.mat', i)));
    rf = rf + S.img;
end
end

function rf = load_foc(P, tag, xf)
rf = zeros(numel(P.z), numel(xf), 'single');
for k = 1:numel(xf)
    S = load(fullfile(P.dataDir,'foc',sprintf('%s_l%03d.mat', tag, k)));
    rf(:,k) = S.img;
end
end

function rf = load_foc_noise(P, xf)
rf = zeros(numel(P.z), numel(xf), 'single');
for k = 1:numel(xf)
    S = load(fullfile(P.dataDir,'noise',sprintf('NF_l%03d.mat', k)));
    rf(:,k) = S.img;
end
end

function [zz, snr] = snr_curve(rf, nrf, P, dx)
ps = double(rf).^2;  pn = double(nrf).^2;
h1 = ones(round(5e-3/P.dz),1); h1 = h1/sum(h1);
h2 = ones(1, round(5e-3/dx));  h2 = h2/sum(h2);
ps = imfilter(imfilter(ps, h1, 'replicate'), h2, 'replicate');
pn = imfilter(imfilter(pn, h1, 'replicate'), h2, 'replicate');
S  = 10*log10(ps ./ max(pn, realmin));
zz = P.z;
snr = median(S, 2);
end

function write_csv(outdir, names, R)
fid = fopen(fullfile(outdir,'table_resolution.csv'),'w');
fprintf(fid, 'method,');
for q = 1:numel(R.zPts), fprintf(fid,'w10_%dmm,w6_%dmm,', round(R.zPts(q)*1e3), round(R.zPts(q)*1e3)); end
fprintf(fid,'\n');
for k = 1:numel(names)
    fprintf(fid,'%s,', names{k});
    for q = 1:numel(R.zPts)
        fprintf(fid,'%.3f,%.3f,', R.w10(k,q)*1e3, R.w6(k,q)*1e3);
    end
    fprintf(fid,'\n');
end
fprintf(fid,'FOCUSED,');
for q = 1:numel(R.zPts), fprintf(fid,'%.3f,%.3f,', R.w10F(q)*1e3, R.w6F(q)*1e3); end
fprintf(fid,'\n');
fclose(fid);

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

fid = fopen(fullfile(outdir,'table_snr.csv'),'w');
fprintf(fid,'method,'); fprintf(fid,'%dmm,', round(R.zdep*1e3)); fprintf(fid,'\n');
for k = 1:numel(names)
    fprintf(fid,'%s,', names{k});
    fprintf(fid,'%.2f,', R.snrRel(k,:));
    fprintf(fid,'\n');
end
fprintf(fid,'FOCUSED,'); fprintf(fid,'%.2f,', R.snrRelF); fprintf(fid,'\n');
fclose(fid);

fid = fopen(fullfile(outdir,'table_speckle.csv'),'w');
fprintf(fid,'method,mu_over_sigma,corrX_mm,corrZ_mm\n');
for k = 1:numel(names)
    fprintf(fid,'%s,%.3f,%.3f,%.3f\n', names{k}, R.speck(k).ratio, ...
        R.speck(k).corrX*1e3, R.speck(k).corrZ*1e3);
end
fprintf(fid,'FOCUSED,%.3f,%.3f,%.3f\n', R.speckF.ratio, R.speckF.corrX*1e3, R.speckF.corrZ*1e3);
fclose(fid);
end
