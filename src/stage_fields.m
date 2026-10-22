function stage_fields(P)
%STAGE_FIELDS  Transmit pressure maps for every plane wave + synthetic field.
%   Uses calc_hp on a 1 mm grid (41 x 56 points).  The synthetic (compounded)
%   transmit field is the coherent sum of the 71 pressure fields, each shifted
%   by its own transmit delay tau_ec(alpha,x,z) + dt(alpha).

fprintf('\n===== transmit field maps =====\n');
C = load(fullfile(P.dataDir, 'calib.mat'));
outdir = fullfile(P.dataDir, 'field');
if ~exist(outdir, 'dir'), mkdir(outdir); end

x = P.fmX;  z = P.fmZ;
[Xg, Zg] = meshgrid(x, z);
npts = numel(Xg);
pts  = [Xg(:), zeros(npts,1), Zg(:)];

Th = f2_setup(P);

K0 = 1200;  K1 = 1200;          % zero padding for the per-pixel time shifts
Nt = [];
syn = [];
peakAll = zeros(numel(z), numel(x), P.nAngles, 'single');
tauRef = 30e-3 / P.c;

for i = 1:P.nAngles
    a = P.angles(i);
    d = P.xe * sin(a) / P.c;
    xdc_center_focus(Th, [0 0 0]);
    xdc_focus_times(Th, 0, d);
    xdc_apodization(Th, 1, ones(1, P.N));

    [hp, t0f] = calc_hp(Th, pts);

    if isempty(Nt)
        Nt  = size(hp, 1);
        syn = zeros(Nt, npts);
    end
    if size(hp,1) < Nt
        hp(end+1:Nt, :) = 0;
    elseif size(hp,1) > Nt
        hp = hp(1:Nt, :);
    end

    peakAll(:,:,i) = single(reshape(max(abs(hp), [], 1), numel(z), numel(x)));

    % Synthetic transmit field: every angle is shifted by its transit time to
    % the FOCUS (0, 30 mm) -- a constant per angle.  Applying instead a
    % per-pixel beamforming shift would align the pulses everywhere and the
    % sum would no longer show a focal spot (it would just be n x the field).
    tauS = 30e-3 * cos(a) / P.c + C.dt(i);
    for c0 = 1:256:npts
        cc = c0:min(npts, c0+255);
        nc = numel(cc);
        % each calc_hp record starts at its own t0f: express the shift in
        % absolute time so the records are aligned across angles
        sh = (tauS - t0f - tauRef) * P.fs * ones(1, nc);

        H  = [zeros(K0,nc); hp(:,cc); zeros(K1,nc)];
        Hr = size(H,1);
        q  = (0:Nt-1).' + sh;                      % Nt x nc
        i0 = floor(q);
        fr = q - i0;
        col = (0:nc-1) * Hr;                       % linear-index column offset
        li = i0 + K0 + 1;
        li = min(max(li, 1), Hr-1) + col;

        vv = H(li).*(1-fr) + H(li+1).*fr;
        syn(:, cc) = syn(:, cc) + vv;
    end

    if mod(i,10) == 0
        fprintf('  field %2d/%d\n', i, P.nAngles);
    end
end

synPeak = reshape(max(abs(syn), [], 1), numel(z), numel(x));

% ---- focused transmit field for comparison (same 52-element aperture) ------
xdc_center_focus(Th, [0 0 0]);
xdc_focus(Th, 0, [0 0 P.txFocus]);
apo = double(abs(P.xe) <= P.txApodW/2 + 1e-12);
xdc_apodization(Th, 1, apo);
[hpF, ~] = calc_hp(Th, pts);
if size(hpF,1) < Nt, hpF(end+1:Nt,:) = 0; else, hpF = hpF(1:Nt,:); end
focPeak = reshape(max(abs(hpF), [], 1), numel(z), numel(x));

% ---- focus gain -----------------------------------------------------------
[~, iz] = min(abs(z - 30e-3));  [~, ix] = min(abs(x - 0));
p0   = squeeze(peakAll(iz, ix, :));
gain = synPeak(iz, ix) / mean(p0);

S = struct('x', x, 'z', z, 'peakAll', peakAll, 'synPeak', synPeak, ...
           'focPeak', focPeak, 'gain', gain, 'Nt', Nt);
save(fullfile(outdir, 'maps.mat'), '-struct', 'S', '-v7.3');
fprintf('saved field maps, focus gain = %.1f (theory %d)\n', gain, P.nAngles);
end
