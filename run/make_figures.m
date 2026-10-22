function make_figures()
%MAKE_FIGURES  Produce F1-F13 from data/calib.mat, data/field/maps.mat and
%   data/res/results.mat.  Called after stage_analysis.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);
C = load(fullfile(P.dataDir, 'calib.mat'));
M = load(fullfile(P.dataDir, 'field', 'maps.mat'));
R = load(fullfile(P.dataDir, 'res', 'results.mat'));

F = P.figDir;
if ~exist(F, 'dir'), mkdir(F); end

% --------------------------------------------------------------------- F1
fig = figure('Position', [100 100 900 480]); hold on;
for k = 1:4:P.N
    rectangle('Position', [(P.xe(k)-P.width/2)*1e3, -0.35, P.width*1e3, 0.35], ...
        'FaceColor', [0.35 0.35 0.35], 'EdgeColor', 'none');
end
xa = linspace(-25, 25, 2);
for a = [-16.5 0 16.5]*pi/180
    x0 = -25*sin(a);
    plot([x0 x0+50*sin(a)], [0 60*cos(a)], 'LineWidth', 1.2);
end
plot([0 0], [0 60], 'k--');
text(-24, 4, 'array, 128 el., 0.33 mm pitch', 'FontSize', 9);
text(6, 55, '\alpha = +16.5\circ', 'FontSize', 9);
text(6, 57, '\alpha = 0\circ', 'FontSize', 9);
text(-24, 55, '\alpha = -16.5\circ', 'FontSize', 9);
xlabel('x (mm)'); ylabel('z (mm)'); axis([-25 25 -2 62]); grid on;
title('Plane-wave transmit geometry');
save_fig(fig, fullfile(F, 'F1_geometry.png')); close(fig);

% --------------------------------------------------------------------- F2
fig = figure('Position', [100 100 900 340]);
subplot(1,2,1);
tw = (0:numel(C.pulseWav)-1)'/P.fs*1e6;
plot(tw, C.pulseWav / max(abs(C.pulseWav)), 'LineWidth', 1.2);
xlabel('t (\mus)'); ylabel('normalised'); grid on;
title(sprintf('Two-way pulse, f_c = %.2f MHz', C.pulse.fc/1e6));
subplot(1,2,2);
plot(C.pulseF/1e6, C.pulseS, 'LineWidth', 1.2); hold on;
yl = ylim; plot([C.pulse.f1 C.pulse.f2]/1e6, [-6 -6], 'r-', 'LineWidth', 1.5);
plot([C.pulse.fc C.pulse.fc]/1e6, [-30 0], 'k:');
xlim([2 8]); ylim([-45 2]); grid on;
xlabel('f (MHz)'); ylabel('dB');
title(sprintf('-6 dB BW = %.1f %%', 100*C.pulse.bw));
save_fig(fig, fullfile(F, 'F2_pulse.png')); close(fig);

% --------------------------------------------------------------------- F3
i0 = find(P.iAng == 0);
mp = double(M.peakAll(:,:,i0)); mp = mp / max(mp(:));
fig = figure('Position', [100 100 950 620]);
subplot(2,2,[1 3]);
imagesc(M.x*1e3, M.z*1e3, 20*log10(max(mp, 1e-4))); axis image; colormap(gca, hot);
caxis([-30 0]); colorbar; xlabel('x (mm)'); ylabel('z (mm)');
title('Transmit pressure, \alpha = 0\circ (dB)');
subplot(2,2,2);
[~, iz] = min(abs(M.z - 30e-3));
plot(M.x*1e3, 20*log10(max(mp(iz,:), 1e-4)), 'LineWidth', 1.2);
xlabel('x (mm)'); ylabel('dB'); grid on; title('lateral at z = 30 mm');
subplot(2,2,4);
[~, ix] = min(abs(M.x));
plot(M.z*1e3, 20*log10(max(mp(:,ix), 1e-4)), 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('dB'); grid on; title('axial at x = 0');
save_fig(fig, fullfile(F, 'F3_field_pw0.png')); close(fig);

% --------------------------------------------------------------------- F4
isel = [1 10 19 28 36 44 53 62 71];
fig = figure('Position', [100 100 900 700]);
for q = 1:9
    subplot(3,3,q);
    mp = double(M.peakAll(:,:,isel(q))); mp = mp / max(mp(:));
    imagesc(M.x*1e3, M.z*1e3, 20*log10(max(mp, 1e-4))); axis image;
    colormap(gca, hot); caxis([-30 0]); set(gca, 'XTick', -20:20:20);
    title(sprintf('\\alpha = %.1f\\circ', P.angles(isel(q))*180/pi), 'FontSize', 8);
    if q == 1, ylabel('z (mm)'); end
    if q >= 7, xlabel('x (mm)'); end
end
save_fig(fig, fullfile(F, 'F4_field_angles.png')); close(fig);

% --------------------------------------------------------------------- F5
syn = double(M.synPeak); syn = syn / max(syn(:));
foc = double(M.focPeak); foc = foc / max(foc(:));
pw1 = double(M.peakAll(:,:,i0)); pw1 = pw1 / max(pw1(:));
fig = figure('Position', [100 100 950 620]);
subplot(2,2,[1 3]);
imagesc(M.x*1e3, M.z*1e3, 20*log10(max(syn, 1e-4))); axis image; colormap(gca, hot);
caxis([-30 0]); colorbar; xlabel('x (mm)'); ylabel('z (mm)');
title(sprintf('Synthetic transmit field, %d angles (focus gain %.0f)', ...
    P.nAngles, M.gain));
subplot(2,2,2);
plot(M.x*1e3, 20*log10(max(syn(iz,:),1e-4)), 'LineWidth', 1.4); hold on;
plot(M.x*1e3, 20*log10(max(foc(iz,:),1e-4)), '--', 'LineWidth', 1.2);
plot(M.x*1e3, 20*log10(max(pw1(iz,:),1e-4)), ':', 'LineWidth', 1.2);
xlabel('x (mm)'); ylabel('dB (own max)'); grid on;
legend('71-angle sum','focused','single PW'); title('lateral at z = 30 mm');
subplot(2,2,4);
plot(M.z*1e3, 20*log10(max(syn(:,ix),1e-4)), 'LineWidth', 1.4); hold on;
plot(M.z*1e3, 20*log10(max(foc(:,ix),1e-4)), '--', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('dB (own max)'); grid on;
legend('71-angle sum','focused'); title('axial at x = 0');
save_fig(fig, fullfile(F, 'F5_field_synthetic.png')); close(fig);

% --------------------------------------------------------------------- F6
S = load(fullfile(P.dataDir, 'pw', 'P3_a036.mat'));
fig = figure('Position', [100 100 1000 420]);
subplot(1,3,[1 2]);
tt = (0:size(S.ch,1)-1)/P.fs*1e6;
imagesc(tt, 1:P.N, 20*log10(abs(double(S.ch)).'/max(abs(double(S.ch(:))))));
axis xy; colormap(gca, gray); caxis([-60 0]);
xlabel('t (\mus)'); ylabel('element'); title('128-channel raw RF (\alpha = 0^\circ)');
subplot(1,3,3);
for e = [16 64 112]
    plot(tt, double(S.ch(:,e))/max(abs(double(S.ch(:)))), 'LineWidth', 0.8); hold on;
end
xlim([30 50]); grid on; xlabel('t (\mus)'); ylabel('normalised');
legend('el 16','el 64','el 112'); title('three channels');
save_fig(fig, fullfile(F, 'F6_channel_rf.png')); close(fig);

% --------------------------------------------------------------------- F7
fig = figure('Position', [100 100 950 340]);
ks = [1 2 7];
for q = 1:3
    subplot(1,3,q);
    rr = double(R.P1.rf(:,:,ks(q)));
    rr = rr / max(abs(rr(:)));
    imagesc(P.x*1e3, P.z*1e3, 20*log10(max(abs(rr),1e-3)));
    axis image; colormap(gca, gray); caxis([-40 0]);
    xlim([-15 15]); ylim([18 45]);
    xlabel('x (mm)'); ylabel('z (mm)');
    title(sprintf('%s (RF before envelope)', R.names{ks(q)}));
end
save_fig(fig, fullfile(F, 'F7_rf_images.png')); close(fig);

% --------------------------------------------------------------------- F8
fig = figure('Position', [100 100 1100 760]);
for q = 1:7
    subplot(2,4,q);
    imagesc(P.x*1e3, P.z*1e3, R.P2.db(:,:,q)); axis image;
    colormap(gca, gray); caxis([-P.drDB 0]);
    title(sprintf('%s (n=%d)', R.names{q}, R.nAngles(q)), 'FontSize', 9);
    if q == 1, ylabel('z (mm)'); end
    if q >= 5, xlabel('x (mm)'); end
end
subplot(2,4,8);
imagesc(R.FOC.x*1e3, P.z*1e3, R.FOC.db2); axis image; colormap(gca, gray);
caxis([-P.drDB 0]); title('focused (128 lines)', 'FontSize', 9); xlabel('x (mm)');
save_fig(fig, fullfile(F, 'F8_bmode_cysts.png')); close(fig);

% --------------------------------------------------------------------- F9
fig = figure('Position', [100 100 1100 760]);
for q = 1:7
    subplot(2,4,q);
    imagesc(P.x*1e3, P.z*1e3, R.P1.db(:,:,q)); axis image;
    colormap(gca, gray); caxis([-P.drDB 0]);
    title(sprintf('%s (n=%d)', R.names{q}, R.nAngles(q)), 'FontSize', 9);
    if q == 1, ylabel('z (mm)'); end
    if q >= 5, xlabel('x (mm)'); end
end
subplot(2,4,8);
imagesc(R.FOC.x*1e3, P.z*1e3, R.FOC.db1); axis image; colormap(gca, gray);
caxis([-P.drDB 0]); title('focused (128 lines)', 'FontSize', 9); xlabel('x (mm)');
save_fig(fig, fullfile(F, 'F9_bmode_uniform.png')); close(fig);

% --------------------------------------------------------------------- F10
fig = figure('Position', [100 100 900 380]);
subplot(1,2,1);
plot(R.zPts*1e3, R.w10.'*1e3, 'o-', 'LineWidth', 1.2); hold on;
plot(R.zPts*1e3, R.w10F*1e3, 'k--s', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('-10 dB lateral width (mm)'); grid on;
legend([R.names, {'focused'}], 'Location', 'northwest', 'FontSize', 7);
title('Lateral resolution');
subplot(1,2,2);
plot(R.zPts*1e3, R.w6.'*1e3, 'o-', 'LineWidth', 1.2); hold on;
plot(R.zPts*1e3, R.w6F*1e3, 'k--s', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('-6 dB axial width (mm)'); grid on;
title('Axial resolution');
save_fig(fig, fullfile(F, 'F10_resolution.png')); close(fig);

% --------------------------------------------------------------------- F11
fig = figure('Position', [100 100 900 380]);
subplot(1,2,1);
Cc = [arrayfun(@(s) s.Cmean, R.contrast), R.contrastF.Cmean];
Ce = [arrayfun(@(s) s.Cstd, R.contrast), R.contrastF.Cstd];
errorbar([R.nAngles 128], Cc, Ce, 'o-', 'LineWidth', 1.2); hold on;
Cn = [arrayfun(@(s) s.Cmean, R.contrastN), R.contrastFN.Cmean];
errorbar([R.nAngles 128], Cn, Cn*0 + [arrayfun(@(s) s.Cstd, R.contrastN), R.contrastFN.Cstd], 's--');
set(gca, 'XScale', 'log'); xlabel('n angles'); ylabel('|C| (dB)'); grid on;
legend('noise-free','with noise','Location','southeast');
title('Cyst contrast');
subplot(1,2,2);
errorbar([R.nAngles 128], [arrayfun(@(s) s.CNRmean, R.contrast), R.contrastF.CNRmean], ...
    [arrayfun(@(s) s.CNRstd, R.contrast), R.contrastF.CNRstd], 'o-', 'LineWidth', 1.2);
set(gca, 'XScale', 'log'); xlabel('n angles'); ylabel('CNR'); grid on;
title('CNR');
save_fig(fig, fullfile(F, 'F11_contrast.png')); close(fig);

% --------------------------------------------------------------------- F12
fig = figure('Position', [100 100 700 420]);
plot(R.zdep*1e3, R.snrRel.', 'o-', 'LineWidth', 1.2); hold on;
plot(R.zdep*1e3, R.snrRelF, 'k--s', 'LineWidth', 1.2);
for q = 1:numel(R.nAngles)
    plot(R.zdep([1 end])*1e3, 10*log10(R.nAngles(q))*[1 1], ':', 'Color', [0.6 0.6 0.6]);
end
xlabel('z (mm)'); ylabel('SNR relative to PW at 30 mm (dB)'); grid on;
legend([R.names, {'focused'}], 'Location', 'northeast', 'FontSize', 7);
title('SNR vs depth (relative)');
save_fig(fig, fullfile(F, 'F12_snr.png')); close(fig);

% --------------------------------------------------------------------- F13
sp = R.speck(1);
fig = figure('Position', [100 100 1000 380]);
subplot(1,3,1);
bw = sp.histEdges(2) - sp.histEdges(1);
bar(sp.histEdges(1:end-1)+bw/2, sp.hist/bw, 1, 'FaceColor', [0.4 0.4 0.4]); hold on;
xx = linspace(0.01, 4, 200);
plot(xx, exp(-xx), 'r-', 'LineWidth', 1.2);
xlabel('I / <I>'); ylabel('pdf'); grid on; legend('measured','exponential (ideal)');
title(sprintf('Intensity histogram, \\mu/\\sigma = %.2f (ideal 1.00)', sp.ratio));
subplot(1,3,2);
for q = 1:7
    plot(0:numel(R.speck(q).hist)-1, R.speck(q).hist, '-', 'LineWidth', 1); hold on;
end
xlabel('bin'); ylabel('fraction per bin'); grid on; title('histograms, all methods');
subplot(1,3,3);
plot(R.nAngles, [R.speck.ratioA], 'o-', 'LineWidth', 1.2); hold on;
plot([1 71], [1.91 1.91], 'k--');
set(gca, 'XScale', 'log'); xlabel('n angles'); ylabel('envelope \mu/\sigma');
grid on; title('speckle statistics vs n');
save_fig(fig, fullfile(F, 'F13_speckle.png')); close(fig);

fprintf('figures written to %s\n', F);
end

function v = mp_lat(M, i0, iz)
v = double(M.peakAll(iz, :, i0));
end
