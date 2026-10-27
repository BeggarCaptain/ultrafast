function c = metrics_contrast(env, x, z, P)
%METRICS_CONTRAST  九个无回声囊肿的对比度与 CNR。
%
%   C = METRICS_CONTRAST(ENV, X, Z, P)
%
%   输入
%     env ：包络图像 [numel(z) x numel(x)]
%     x,z ：坐标向量 [m]
%     P   ：参数结构体（提供囊肿中心与半径）
%
%   输出结构体（每个囊肿一项，另有九个囊肿的均值/标准差）
%     C, CNR        ：逐囊肿的对比度 [dB] 与对比度噪声比
%     Cmean, Cstd   ：|C| 在九个囊肿上的均值与标准差 [dB]
%     CNRmean, CNRstd
%
%   定义（设计文档式 6）
%     C   = 10 log10(mu_in / mu_out)   mu 取平均“强度”（包络的平方）
%     CNR = |mu_out - mu_in| / sqrt(var_out + var_in)
%
%   为什么是 10log10 而不是 20log10：mu 是强度（功率），而 dB 的通用刻度是
%   10log10(功率比)；它与 20log10(幅度比) 完全等价（因为 20log10 A = 10log10 A^2），
%   也正是 to_db 里对包络取 20log10 的那条刻度。若对强度比直接取 20log10，
%   得到的 dB 值会正好是正确值的两倍——这是本项目曾经犯过的一个错误。
%
%   区域：每个囊肿内部半径 2 mm 的圆盘（完全无回声），外部 4–6 mm 的环带
%   （背景散斑）。两区之间留 2 mm 间隔，是为了把囊壁和部分容积效应排除在
%   统计之外。
%
%   为什么用强度而不是包络：论文的对比度定义在强度图上；而且散斑的强度服从
%   指数分布，CNR 的归一化正是基于这个假设。
%
%   注意（报告 5.6）：囊内接近全黑时，CNR 会饱和在 0.7 附近，不再能区分方法
%   优劣；在本文这种“无噪声/旁瓣受限”的条件下，有判别力的是 |C|。

I = double(env).^2;                     % 强度 = 包络的平方
[XX, ZZ] = meshgrid(x, z);
Cx = cyst_centres(P);                   % 九个中心（3 列 x 3 行）
n  = size(Cx, 1);

C = zeros(1,n);  CNR = zeros(1,n);
for k = 1:n
    r  = sqrt((XX - Cx(k,1)).^2 + (ZZ - Cx(k,2)).^2);   % 到该囊肿中心的距离
    in = r <= 2e-3;                                     % 内部圆盘：无回声
    ou = r >= 4e-3 & r <= 6e-3;                         % 外部环带：散斑

    mi = mean(I(in));  mo = mean(I(ou));                % 内部/外部平均强度
    vi = var(I(in));   vo = var(I(ou));                 % 内部/外部方差

    C(k)   = 10*log10(mi/mo);                           % 对比度，dB（负值）
    CNR(k) = abs(mo - mi) / sqrt(vo + vi);              % 对比度噪声比
end

c.C = C;  c.CNR = CNR;
c.Cmean = mean(abs(C));  c.Cstd = std(abs(C));          % 取绝对值，报告为正值
c.CNRmean = mean(CNR);   c.CNRstd = std(CNR);
end
