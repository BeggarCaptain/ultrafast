function acquire_all_pw(P, tag, pos, amp, keepCh)
%ACQUIRE_ALL_PW  Acquire + beamform all 71 plane-wave angles for one phantom.
%   Files: data/pw/<TAG>_a###.mat  (resume-safe: existing files are skipped)

outdir = fullfile(P.dataDir, 'pw');
if ~exist(outdir, 'dir'), mkdir(outdir); end
C  = load(fullfile(P.dataDir, 'calib.mat'));
Th = f2_setup(P);

tAll = tic;
for i = 1:P.nAngles
    f = fullfile(outdir, sprintf('%s_a%03d.mat', tag, i));
    if exist(f, 'file'), continue; end

    t1 = tic;
    [rf, t0] = acquire_pw(Th, P, P.angles(i), pos, amp);
    dt = C.dt(i);
    [img, w2] = das_pw(rf, t0, P.angles(i), P, P.x, P.z, dt);

    S = struct('img', single(img), 'w2', single(w2), 't0', t0, ...
               'dt', dt, 'alpha', P.angles(i), 'tag', tag);
    if keepCh || (strcmp(tag,'P3') && i == 36)   % centre angle for F6
        S.ch = single(rf);
    end
    save(f, '-struct', 'S');

    el = toc(tAll);
    fprintf('[%s] angle %2d/%d  %5.1f s  elapsed %5.1f min  eta %5.1f min\n', ...
        tag, i, P.nAngles, toc(t1), el/60, el/60*(P.nAngles-i)/i);
end
end
