function run_test_p3()
%RUN_TEST_P3  Quick end-to-end check on the point-target phantom (P3).
%   Validates E0 + DAS + focused transmit before the multi-hour P1 run.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'src'));
P = params(root);
if ~exist(fullfile(P.dataDir,'logs'),'dir'), mkdir(fullfile(P.dataDir,'logs')); end
diary(fullfile(P.dataDir,'logs','test_p3.log')); diary on;
fprintf('=== P3 test %s ===\n', datestr(now));

if ~exist(fullfile(P.dataDir,'calib.mat'),'file')
    stage0_calib(P);
end

[pp, pa] = phantom_points(P);
acquire_all_pw(P, 'P3', pp, pa, false);
acquire_all_focused(P, 'P3', pp, pa, P.xpsf, false);

% ---- assemble and measure -------------------------------------------------
Sp = load(fullfile(P.dataDir,'pw','P3_a036.mat'));
env = abs(hilbert(double(Sp.img)));
m = metrics_resolution(env, P.x, P.z, 0, 30e-3);
fprintf('PW 0 deg : peak (%.2f, %.2f) mm, w10 %.2f mm, w6 %.2f mm\n', ...
    m.xpk*1e3, m.zpk*1e3, m.w10*1e3, m.w6*1e3);

foc = zeros(numel(P.z), numel(P.xpsf));
for k = 1:numel(P.xpsf)
    S = load(fullfile(P.dataDir,'foc',sprintf('P3_l%03d.mat',k)));
    foc(:,k) = S.img;
end
envF = abs(hilbert(double(foc)));
mF = metrics_resolution(envF, P.xpsf, P.z, 0, 30e-3);
fprintf('focused  : peak (%.2f, %.2f) mm, w10 %.2f mm, w6 %.2f mm\n', ...
    mF.xpk*1e3, mF.zpk*1e3, mF.w10*1e3, mF.w6*1e3);

% ---- apodisation check ----------------------------------------------------
Sf = load(fullfile(P.dataDir,'foc','P3_l128.mat'));
if isfield(Sf,'ch')
    e = mean(abs(Sf.ch),1);
    fprintf('transmit apodisation: edge/centre channel energy ratio %.3f\n', ...
        mean(e([1:3 end-2:end]))/mean(e(60:68)));
end

% ---- figures --------------------------------------------------------------
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
