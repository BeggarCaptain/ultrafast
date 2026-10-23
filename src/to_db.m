function B = to_db(A, dr)
%TO_DB  把包络（或任意非负图）转换成以全局最大值为参考的 dB 图，并从下方截断。
%
%   B = TO_DB(A, DR) 返回 20*log10(A / max(A))，数值小于 -DR 的部分被截到 -DR。
%   参考值取 A 的全局最大值，因此若多幅图要共用同一色标，请调用前自行统一归一
%   化（本项目的复合图都是按各自最大值归一化后再比较形状的）。
%
%   输入
%     A  ：非负数组（Hilbert 包络、|RF|、声压幅值等）
%     DR ：显示动态范围，单位 dB（例如 60）；传 [] 或省略表示不截断
%
%   输出
%     B  ：同尺寸数组，取值范围 [-DR, 0]
%
%   用 realmin 兜底是为了避免 log10(0) = -Inf：一旦出现 -Inf，后续任何运算都
%   会变成 NaN，图像色标就会整幅失效。

A = double(A);                              % single 转 double：取对数需要足够动态范围
B = 20*log10(max(A, realmin) / max(A(:)));  % 先按全局峰值归一化，再取 20log10
if nargin > 1 && ~isempty(dr)
    B = max(B, -dr);                        % 噪声底截断到显示动态范围
end
end
