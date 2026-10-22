function acquire_all_focused(P, tag, pos, amp, xf, keepCh)
%ACQUIRE_ALL_FOCUSED  Acquire + beamform every focused line for one phantom.
%   Files: data/foc/<TAG>_l###.mat
%   XF : vector of line positions [m]

outdir = fullfile(P.dataDir, 'foc');
if ~exist(outdir, 'dir'), mkdir(outdir); end
C  = load(fullfile(P.dataDir, 'calib.mat'));
Th = f2_setup(P);

nL = numel(xf);
tAll = tic;
for k = 1:nL
    f = fullfile(outdir, sprintf('%s_l%03d.mat', tag, k));
    if exist(f, 'file'), continue; end

    t1 = tic;
    [rf, t0] = acquire_focused(Th, P, xf(k), pos, amp);
    [img, w2] = das_focused(rf, t0, xf(k), P, P.z, C.dt_foc);

    S = struct('img', single(img), 'w2', single(w2), 't0', t0, ...
               'dt', C.dt_foc, 'xf', xf(k), 'tag', tag, 'Nt', size(rf,1));
    if keepCh && k == round(nL/2)
        S.ch = single(rf);
    end
    save(f, '-struct', 'S');

    el = toc(tAll);
    if mod(k, 10) == 0 || k == nL
        fprintf('[%s] line %3d/%d  %5.1f s  elapsed %5.1f min  eta %5.1f min\n', ...
            tag, k, nL, toc(t1), el/60, el/60*(nL-k)/k);
    end
end
end
