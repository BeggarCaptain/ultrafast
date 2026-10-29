function stage0_calib(P)
%STAGE0_CALIB  E0：脉冲、延时模型、时间零点、孔径、插值等自检。
%
%   STAGE0_CALIB(P) 在正式长采集之前完成必须的自检，并把后续各阶段需要的
%   东西写入 data/calib.mat（每个角度的 dt、聚焦模式的 dt、脉冲指标等）。
%
%   为什么需要标定阶段？
%   --------------------
%   Field II 与自写 DAS 必须对“时间零点”有一致的理解。通道数据的时间轴由
%   Field II 给出（含换能器冲激响应与滤波带来的群延迟），DAS 用的是自己算的
%   几何往返时间；两者相差一个常数。常数偏移只让整幅图在深度上平移（
%   dz = c*dt/2），但**每个角度各自的偏移不同**时，图像会被剪切
%   （dx = -z*sin(alpha)），相干复合直接失效。下面的检查是把这两种效应量出来，
%   而不是假设它们不存在。
%
%   dt 是怎么来的、为什么乘 2
%   ------------------------
%   把一个点散射子放在 30 mm，先令 dt = 0 成像，峰值会落在 z_obs。模型延时比
%   真实往返时间少了（或多了）2*(z_obs - 30 mm)/c：延时误差在去、回两程各出现
%   一次，所以要把深度差乘 2 再除声速。把 dt 加进模型延时后，峰值回到 30 mm。
%   dt 对所有像素是同一个数——这正是它可以用一个标量表示的原因。
%
%   内容
%     (1) 脉冲频谱     -> 两程 fc、-6 dB 带宽、半功率带宽（约 4.1 MHz；半功率约 62 %）
%     (2) 平面波倾角   -> 峰值时刻随 x 的斜率（应为 sin(a)/c）
%     (3) 时间零点     -> 每个角度一个 dt，聚焦模式另有一个
%     (4) 点目标 PSF   -> 30 mm 处的横向 / 轴向 -6 dB 宽度
%     (5) 接收孔径     -> 各深度的阵元数（F 数 = 1.75）
%     (6) 线性/样条插值 -> 插值误差有多大
%     (7) 采集耗时     -> 一次 2 万散射子发射要多久
%
%   验收标准（报告里写“通过”的依据；括号内是本机实测值）
%     发射载频    4.50 MHz                   参数给定（报告 3.1 节），不是被测对象
%     两程 fc     约 4.1 MHz                 点散射子实测（4.12 MHz）；孔径使其比
%                                            发射低约 8 %，只记录、不校正
%     半功率带宽  60 ± 5 %                   同一测量（62.5 %）；60 % 是按半功率调的
%     -6 dB 带宽  仅报告，无目标             同一频谱的另一个电平（86 % 量级）
%     tilt ratio  1 ± 0.05（5° 与 10°）      0° 斜率为零，不参与比值（实测 1.000/1.000）
%     dt 极差     ≤ 0.1 µs                   更大说明模型有系统性错误（实测 0.014 µs）
%     峰值深度    校正后 30.00 ± 0.05 mm，四个测试角度（实测 30.00 mm）
%     PSF         30 mm 处 wLat ≈ 0.85 ± 0.1 mm，wAx ≈ 0.2 ± 0.05 mm
%                 （实测 0.83 / 0.19 mm；注意这是**单平面波**，0.7 mm 是 11 角度
%                  复合后的值，两者不能混用）
%     接收孔径    10/30/60 mm 分别为 18/52/104 阵元
%     插值误差    相对峰值低于 -20 dB（实测 -22.6 dB）
%
%   任何一项不通过，都不要开始长采集。常见原因：延时符号写反、DAS 里漏掉
%   +x1 sin(alpha)/c 项、脉冲形状没调到目标带宽。
%
%   输出：data/calib.mat（结构体 S，见函数末尾）。后续所有采集都从这里读 dt，
%   所以标定只需做一次；改动采集相关参数后它会被 stamp_guard 判为过期。

fprintf('\n===== E0 calibration =====\n');
[Th, ThR] = f2_setup(P);                       % 发射句柄 + 接收句柄
ptpos = [0 0 30e-3];  ptamp = 1;             % 30 mm 处一个弱散射子

% ---------------------------------------------------------------- (1) 脉冲
% 轴线上的点散射子，每个通道拿到的就是两程脉冲。
[rf1, t01] = acquire_pw(Th, ThR, P, 0, ptpos, ptamp);
fprintf('single-scatterer record: %d samples x %d channels, t0 = %.3f us\n', ...
    size(rf1,1), size(rf1,2), t01*1e6);

ps = pulse_stats(rf1, P.sampleRate);
% 这是**两程**脉冲：点散射子的回波在发射、接收各经过一次换能器响应，频域上
% 等于单程谱的平方（再乘孔径的空间响应）。所以它的中心频率比发射低（4.12 对
% 4.50 MHz）、带宽也窄（86 % 对 112 %）。图 F2 画的是单程的理想发射脉冲，
% 两者不能混为一谈。
% 两个电平都打印：-6 dB 与分辨率指标同刻度，半功率则是“60 % 设计目标”的口径
fprintf('pulse: fc = %.3f MHz, -6 dB BW = %.1f %% (half power %.1f %%)\n', ...
    ps.fc/1e6, 100*ps.bw, 100*ps.bwHp);

% ------------------------------------------------------------ (2) 倾角检查
% 40 mm 深度排一行点：若平面波倾斜正确，峰值到达时刻应随 x 线性增加，斜率
% 为 sin(alpha)/c。这一步能抓住“符号写反”或“延时被截断”的问题。
% 为什么取 40 mm、±15 mm：40 mm 已远离近场（阵列半长 21 mm），波前已充分倾斜；
% ±15 mm 覆盖了大部分有效孔径，斜率拟合的信噪比最好。
xchk = (-15:1:15)*1e-3;
pchk = [xchk(:), zeros(numel(xchk),1), 40e-3*ones(numel(xchk),1)];
tilt = zeros(1,3);  slope = zeros(1,3);
alphas_chk = [0 5 10]*pi/180;
for q = 1:3
    d = P.elementX*sin(alphas_chk(q))/P.soundSpeed;
    xdc_center_focus(Th,[0 0 0]); xdc_focus_times(Th, 0, d);
    xdc_apodization(Th, 0, ones(1,P.nElements));     % 时间线 0 s 起生效（不是通道编号）
    [hp, th0] = calc_hp(Th, pchk);
    tt = th0 + (0:size(hp,1)-1)'/P.sampleRate;
    [~, im] = max(abs(hp), [], 1);           % 每个网格点的峰值时刻
    tp = tt(im);
    cf = polyfit(xchk(:), tp(:), 1);         % 线性拟合：斜率 = sin(alpha)/c
    slope(q) = cf(1);
    if abs(sin(alphas_chk(q))) > 1e-9
        tilt(q) = cf(1) / (sin(alphas_chk(q))/P.soundSpeed);   % 实测/理论比值
    else
        tilt(q) = NaN;                       % 0°：斜率为零，比值无定义
    end
    fprintf('alpha=%4.1f deg: measured slope %.4f us/mm, theory %.4f us/mm, ratio %.3f\n', ...
        alphas_chk(q)*180/pi, cf(1)*1e6, sin(alphas_chk(q))/P.soundSpeed*1e6, tilt(q));
end

% --------------------------------------- (3) 逐角度时间零点标定
% dt = 0 时，30 mm 处的点目标会出现在 z_obs；把模型延时整体加上
% delta = 2 (z_obs - 30 mm)/c，峰值就回到真实位置（设计文档 4.2 节的时间零点）。
% 乘 2 是因为延时误差在去、回两程各出现一次。
% 每个角度只做一条线的 DAS：足够便宜，又能把与角度有关的部分从固定群延迟里
% 分离出来。zcal 只取 26–34 mm 这一窄段——目的只是找峰值位置，窗开得越大越
% 容易把旁瓣框进来。
zcal = (26:0.05:34)'*1e-3;
dt = zeros(1, P.nAngles);
for i = 1:P.nAngles
    [rf, t0] = acquire_pw(Th, ThR, P, P.angles(i), ptpos, ptamp);
    img = das_pw(rf, t0, P.angles(i), P, 0, zcal, 0);      % 这里刻意用 dt = 0
    dt(i) = 2*(peak_pos(zcal, abs(hilbert(img))) - 30e-3)/P.soundSpeed;
end
% 剔除孤立的“跳瓣”异常点：dt 本该是 alpha 的平滑函数。peak_pos 已经很稳，
% 但偶发跳一个轴向旁瓣仍会让该角度错位约 1 mm，把复合图像拖花。判据是"偏离
% 左右各两个角度的中位数超过 0.5 µs"——相邻角度只差 2°，dt 不可能真的跳这么多。
for i = 1:P.nAngles
    lo = max(1, i-2);  hi = min(P.nAngles, i+2);
    med = median(dt(lo:hi));
    if abs(dt(i) - med) > 0.5e-6
        fprintf('  dt outlier at angle %d: %.3f -> %.3f us\n', i, dt(i)*1e6, med*1e6);
        dt(i) = med;
    end
end
fprintf('dt(alpha): min %.3f us, max %.3f us, mean %.3f us, spread %.3f us\n', ...
    min(dt)*1e6, max(dt)*1e6, mean(dt)*1e6, (max(dt)-min(dt))*1e6);

% 聚焦参照：传统 B 模式用的是另一种发射延时，时间零点也不同
[rfF, t0F] = acquire_focused(Th, ThR, P, 0, ptpos, ptamp);
imgF = das_focused(rfF, t0F, 0, P, zcal, 0);
dt_foc = 2*(peak_pos(zcal, abs(hilbert(imgF))) - 30e-3)/P.soundSpeed;
fprintf('dt(focused, xf=0) = %.3f us\n', dt_foc*1e6);

% 验证：校正后峰值必须落在 30.00 mm。取三个分布开的角度就足以发现系统性错误。
zver = zeros(1,3);  aver = [1 6 11];
for q = 1:3
    i = aver(q);
    [rf, t0] = acquire_pw(Th, ThR, P, P.angles(i), ptpos, ptamp);
    img = das_pw(rf, t0, P.angles(i), P, 0, zcal, dt(i));
    zver(q) = peak_pos(zcal, abs(hilbert(img)));
end
fprintf('calibrated peak depth (angles 1/6/11): %s mm\n', ...
    mat2str(round(zver*1e3, 3)));

% ------------------------------------------- (4) 校正后的点目标 PSF
% 目标附近一小块二维图像：给出 30 mm 处的横向、轴向 -6 dB 宽度参考值，
% 供 E4 与其他深度比较。
xg = (-5:0.05:5)*1e-3;  zg = (26:0.05:34)'*1e-3;   % 注意是列向量！
[rf, t0] = acquire_pw(Th, ThR, P, 0, ptpos, ptamp);
img = das_pw(rf, t0, 0, P, xg, zg, dt(P.angleIndex==0));
env = to_env(img);
[~, imax] = max(env(:));  [iz, ix] = ind2sub(size(env), imax);
lat = env(iz, :)';  axp = env(:, ix);
fprintf('PW alpha=0 point: peak at (%.2f, %.2f) mm\n', xg(ix)*1e3, zg(iz)*1e3);
fprintf('  -6 dB lateral width %.2f mm, -6 dB axial width %.2f mm\n', ...
    crossing_width(xg, lat, 0.5)*1e3, crossing_width(zg, axp, 0.5)*1e3);

% 聚焦线的横向 PSF 需要多条线，留到 E4 再测。

% -------------------------------------------------------- (5) 孔径阵元数
% F 数 = 1.75 意味着接收孔径随深度增大；DAS 的孔径掩膜必须复现这个阵元数。
% 查 10/30/60 mm 三点就够了：浅端（18 个阵元）、焦点（52 个）、深端（104 个），
% 正是孔径从"被物理孔径限制"过渡到"被 F 数限制"的三段。
for zc = [10 30 60]*1e-3
    fprintf('aperture at z=%2.0f mm: %3d elements (%.1f mm)\n', zc*1e3, ...
        sum(abs(P.elementX) <= zc/(2*P.fNumber)), 2*zc/(2*P.fNumber)*1e3);
end

% ------------------------------------------------------ (6) 线性 vs 样条
% 全流程用线性插值以换取速度；这里量化它相对三次样条在一条 A 线上的代价。
% 为什么值得量：DAS 里每个像素、每个通道都要插值一次，样条比线性慢约 5–10 倍；
% 只有确认误差远低于主瓣（这里 −22.6 dB）才敢一直用线性。
% 判据取"相对峰值的最大偏差"而不是绝对误差：绝对误差随图像增益变，没有意义。
[rf, t0] = acquire_pw(Th, ThR, P, 0, ptpos, ptamp);
zc = P.imageZ;
imgLin = das_pw(rf, t0, 0, P, 0, zc, dt(P.angleIndex==0));
imgSp  = zeros(numel(zc),1);
tt = t0 + (0:size(rf,1)-1)'/P.sampleRate;
for c = 1:size(rf,2)                        % 逐通道做样条插值
    tau = zc/P.soundSpeed + sqrt(zc.^2 + P.elementX(c)^2)/P.soundSpeed + dt(P.angleIndex==0);
    m = abs(P.elementX(c)) <= zc/(2*P.fNumber);
    v = interp1(tt, rf(:,c), tau, 'spline', 0);
    imgSp = imgSp + m.*v;
end
d1 = abs(hilbert(imgLin)); d2 = abs(hilbert(imgSp));
band = zc > 20e-3 & zc < 50e-3;             % 只在有效深度段比较
err = 20*log10(max(abs(d1(band)-d2(band))) / max(d1(band)));
fprintf('linear vs spline beamforming: max difference %.1f dB (20-50 mm)\n', err);

% ----------------------------------------------------------- (7) 耗时测试
% 一次全量发射，用来估计整个采集要跑多久。散射子数与 P1 相同（20000 个），
% 所以这里测出的 tAcq20k 可以直接外推：11 个角度约 5 min、128 条聚焦线约 58 min。
% 这个数进报告的性能表，也是"为什么要把采集拆成断点续跑"的依据。
[pos, amp] = phantom_uniform(P);
tic; [rf20, ~] = acquire_pw(Th, ThR, P, 0, pos, amp); tAcq = toc;
fprintf('calc_scat_multi, %d scatterers, %d ch: %.1f s -> record %d samples\n', ...
    numel(amp), P.nElements, tAcq, size(rf20,1));

% --------------------------------------------------------------- 保存结果
% 后续阶段需要的东西集中到一个结构体。脉冲频谱与波形也留着，这样画 F2 时不必
% 重新跑 Field II。
S = struct();
S.matlabVersion = version;
S.computer      = computer;
S.fieldIIPath   = which('field_init');
S.pulse         = rmfield(ps, {'f','S','wav'});   % 标量脉冲指标
S.pulseF        = ps.f;  S.pulseS = ps.S;         % 两程谱（核对用；图 F2 画的是理想发射脉冲）
S.pulseWav      = ps.wav;                         % 两程波形（同上）
S.t01           = t01;
S.tiltRatio     = tilt;
S.tiltSlope     = slope;
S.alphaChk      = alphas_chk;
S.dt            = dt;                             % 逐角度时间零点
S.dt_foc        = dt_foc;                         % 聚焦模式时间零点
S.pwPeak        = [xg(ix) zg(iz)];
S.pwWidthLat    = crossing_width(xg, lat, 0.5);
S.pwWidthAx     = crossing_width(zg, axp, 0.5);
S.interpErrDB   = err;
S.tAcq20k       = tAcq;
S.recordSamples = size(rf20,1);
S.zc            = zc;
S.imgLin        = imgLin;
S.imgSpline     = imgSp;

if ~exist(P.dataDir,'dir'), mkdir(P.dataDir); end
save(fullfile(P.dataDir,'calib.mat'), '-struct', 'S');
fprintf('saved %s\n', fullfile(P.dataDir,'calib.mat'));
end

function e = to_env(x)
%TO_ENV  Hilbert 包络（本地小函数，避免再开一个文件）。
e = abs(hilbert(x));
end
