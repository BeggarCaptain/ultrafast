function img = das_pw_gpu(rf, t0, alpha, P, dt)
%DAS_PW_GPU  DAS_PW 的 GPU 版本（只被可选的 E7 基准测试使用）。
%
%   IMG = DAS_PW_GPU(RF, T0, ALPHA, P, DT)
%
%   RF 必须已经是 gpuArray。DT 与 das_pw 相同，是该角度的时间零点校正 [s]。
%   本函数是 das_pw.m 的逐行移植：延时模型、孔径
%   掩膜、线性插值完全相同，这样 CPU/GPU 的对比测的是“设备”而不是“另一个
%   算法”。只有两处刻意不同：
%     * 不输出噪声权重 w2（E7 不需要）；
%     * 网格直接取 P（基准测试总是用整幅图像）。
%
%   为什么单独写一个函数，而不是在 das_pw 里加 if isa(rf,'gpuArray')
%   ----------------------------------------------------------------
%   das_pw 是全部正式结果（图像、指标、报告）的生产代码，必须保持简单、可读、
%   在无 GPU 的机器上也能跑。把 GPU 分支塞进去，等于让每一次波束形成都多一层
%   与结果无关的判断，也让"参考实现"和"被测实现"混在一起。分开之后，两边代码
%   可以逐行对照，E7 的结论才有说服力。
%
%   为什么仍然按列循环（没有把整幅图一次性向量化）
%   --------------------------------------------
%   一次算出 Nz x Nx x Nch = 1101 x 256 x 128 的中间量约 144 MB（single），
%   在显存里是一次不必要的高峰；按列循环时中间量只有 Nz x Nch = 0.56 MB，
%   而 GPU 的并行度已经由 1101 x 128 这一维提供。
%
%   关于速度：延时表在每次调用时重新计算，这才是“compute-only”的诚实计时；
%   把延时表预先算好常驻显存是合理的进一步优化，但那样报出的就不是本文的算法
%   速度了。

z = P.imageZ(:);  x = P.imageX(:).';              % 深度为列，横向为行
[Nt, Nch] = size(rf);
Nz = numel(z);  Nx = numel(x);

img = gpuArray.zeros(Nz, Nx, 'single');
ca = cos(alpha);  sa = sin(alpha);
off = gpuArray((0:Nch-1) * Nt);         % 通道偏移，常驻显存
dxe = gpuArray(P.elementX);                   % 阵元坐标，同样常驻显存

for k = 1:Nx                            % 与 CPU 一样，一次一条线
    xc  = x(k);
    dxc = dxe - xc;                     % 阵元到像素的距离（gpuArray）
    hw  = z / (2*P.fNumber);               % 动态接收孔径
    m   = abs(dxc) <= hw;               % Nz x N 孔径掩膜

    % 与 das_pw 完全相同的四项，包含 Field II 发射延时修正项。
    tau = (z*ca + xc*sa)/P.soundSpeed + sqrt(z.^2 + dxc.^2)/P.soundSpeed + dxc*sa/P.soundSpeed + dt;
    pos = (tau - t0) * P.sampleRate + 1;        % 小数采样下标

    i0 = floor(pos);
    fr = pos - i0;
    m  = m & (i0 >= 1) & (i0 < Nt);     % 丢掉超出记录范围的延时
    i0 = min(max(i0, 1), Nt-1);         % 取值前钳位

    v0 = rf(i0 + off);                  % 左侧采样，一次取全部通道
    v1 = rf(i0 + 1 + off);              % 右侧采样
    img(:, k) = sum(((1-fr).*m).*v0 + (fr.*m).*v1, 2);
end
end
