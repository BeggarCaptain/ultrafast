function stage_fields(P)
%STAGE_FIELDS  各平面波的发射声压图，以及“合成发射声场”。
%
%   STAGE_FIELDS(P) 用 calc_hp 在 1 mm 网格（41 x 56 点）上计算：
%     * 每个角度的发射峰值声压                 -> peakAll(:,:,i)
%     * 11 个角度声场的相干叠加                 -> synPeak
%     * 一次真实聚焦发射的声场（作为参照）      -> focPeak
%   结果存入 data/field/maps.mat。
%
%   这里的“合成发射声场”指什么
%   --------------------------
%   聚焦发射会让各阵元的贡献同时到达焦点。相干复合在事后处理中达到同样的效果：
%   把每个角度的声场按其“到焦点的传播时间”平移后相加。焦点处各脉冲对齐（增益
%   约 n 倍），离开焦点就不再对齐，于是出现一个焦斑。图 F4/F5 画的就是这个场。
%
%   重要：平移量是“每个角度一个常数”（到焦点的传播时间），不是逐像素的波束
%   形成延时。若用逐像素延时，脉冲会在每个像素都对齐，叠加结果只是单角度声场
%   的 n 倍，根本看不到焦斑。
%
%   输出用 single 精度；原始声压记录不保存（每个角度约 44 MB）。
%
%   为什么用 calc_hp 而不是 calc_scat_multi
%   ---------------------------------------
%   calc_hp 直接给出指定点上的发射声压波形，用来画“场”；calc_scat_multi 返回的
%   是加了散射子、接收延时与脉冲回波之后的通道信号，反推场反而更绕。两者共用
%   同一套阵列与激励，所以 calc_hp 的结果就是这个系统真实的发射场。
%   代价是每个网格点都要算一次空间冲激响应：2296 个点、一个角度约 0.3 s。

fprintf('\n===== transmit field maps =====\n');
C = load(fullfile(P.dataDir, 'calib.mat'));
outdir = fullfile(P.dataDir, 'field');
if ~exist(outdir, 'dir'), mkdir(outdir); end

x = P.fmX;  z = P.fmZ;                       % 1 mm 网格
[Xg, Zg] = meshgrid(x, z);
npts = numel(Xg);
pts  = [Xg(:), zeros(npts,1), Zg(:)];        % calc_hp 要求 N x 3 的点坐标

Th = f2_setup(P);

K0 = 1200;  K1 = 1200;                       % 时间平移用的零填充长度
Nt = [];                                     % 记录长度，由第一次调用确定
syn = [];                                    % 累加中的相干和
peakAll = zeros(numel(z), numel(x), P.nAngles, 'single');
tauRef = 30e-3 / P.c;                        % 任意常数，只有差值有意义

for i = 1:P.nAngles
    a = P.angles(i);
    d = P.xe * sin(a) / P.c;                 % 与采集时完全相同的发射延时
    xdc_center_focus(Th, [0 0 0]);
    xdc_focus_times(Th, 0, d);
    xdc_apodization(Th, 1, ones(1, P.N));

    [hp, t0f] = calc_hp(Th, pts);            % [Nt x npts] 声压，各自的时间基准

    % ---- 保证每个角度的记录长度一致 ------------------------------------
    if isempty(Nt)
        Nt  = size(hp, 1);
        syn = zeros(Nt, npts);
    end
    if size(hp,1) < Nt
        hp(end+1:Nt, :) = 0;                 % 短的补零
    elseif size(hp,1) > Nt
        hp = hp(1:Nt, :);                    % 长的截断
    end

    % ---- 单角度峰值声压图（供 F3/F4） ----------------------------------
    peakAll(:,:,i) = single(reshape(max(abs(hp), [], 1), numel(z), numel(x)));

    % 合成发射场：每个角度按其到焦点 (0, 30 mm) 的传播时间平移，每个角度一个
    % 常数。若改成逐像素的波束形成延时，脉冲处处对齐，就只剩 n 倍声场、没有焦斑。
    tauS = 30e-3 * cos(a) / P.c + C.dt(i);
    for c0 = 1:256:npts                     % 分块，控制索引矩阵的大小
        cc = c0:min(npts, c0+255);
        nc = numel(cc);
        % 每次 calc_hp 的记录都从各自的 t0f 开始：把平移量换算到绝对时间上，
        % 各角度的记录才能对齐。若不减 t0f，记录会差最多约 0.8 µs，叠加不相干。
        sh = (tauS - t0f - tauRef) * P.fs * ones(1, nc);

        H  = [zeros(K0,nc); hp(:,cc); zeros(K1,nc)];   % 前后补零的列
        Hr = size(H,1);
        q  = (0:Nt-1).' + sh;                      % Nt x nc 的小数下标
        i0 = floor(q);
        fr = q - i0;
        col = (0:nc-1) * Hr;                       % 线性索引的列偏移
        li = i0 + K0 + 1;                          % 注意：必须加 col，
        li = min(max(li, 1), Hr-1) + col;          % 否则每一列都读到第 1 列

        vv = H(li).*(1-fr) + H(li+1).*fr;          % 时间方向的线性插值
        syn(:, cc) = syn(:, cc) + vv;              % 相干累加
    end

    if mod(i,10) == 0
        fprintf('  field %2d/%d\n', i, P.nAngles);
    end
end

synPeak = reshape(max(abs(syn), [], 1), numel(z), numel(x));   % 叠加后的峰值

% ---- 作为对照的聚焦发射声场（同样是 52 阵元孔径） ----------------------
xdc_center_focus(Th, [0 0 0]);
xdc_focus(Th, 0, [0 0 P.txFocus]);
apo = double(abs(P.xe) <= P.txApodW/2 + 1e-12);
xdc_apodization(Th, 1, apo);
[hpF, ~] = calc_hp(Th, pts);
if size(hpF,1) < Nt, hpF(end+1:Nt,:) = 0; else, hpF = hpF(1:Nt,:); end
focPeak = reshape(max(abs(hpF), [], 1), numel(z), numel(x));

% ---- 焦点增益 -----------------------------------------------------------
% 焦点处相干叠加的峰值 ÷ 同一点处单角度峰值的平均值。理论值 n（=11）；实测会
% 略小于 n，因为各角度的脉冲波形并不完全相同（见报告 5.3）。
[~, iz] = min(abs(z - 30e-3));  [~, ix] = min(abs(x - 0));
p0   = squeeze(peakAll(iz, ix, :));
gain = synPeak(iz, ix) / mean(p0);

S = struct('x', x, 'z', z, 'peakAll', peakAll, 'synPeak', synPeak, ...
           'focPeak', focPeak, 'gain', gain, 'Nt', Nt);
save(fullfile(outdir, 'maps.mat'), '-struct', 'S', '-v7.3');   % 逐角度声场，用 v7.3
fprintf('saved field maps, focus gain = %.1f (theory %d)\n', gain, P.nAngles);
end
