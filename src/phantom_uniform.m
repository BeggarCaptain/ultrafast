function [pos, amp] = phantom_uniform(P)
%PHANTOM_UNIFORM  P1：x-z 平面内均匀随机分布的散射子（y = 0）。
%
%   [POS, AMP] = PHANTOM_UNIFORM(P) 返回
%     POS ：P.nScatterers x 3 矩阵 [x y z]，单位米，y = 0（二维仿真）
%     AMP ：P.nScatterers x 1 幅度，零均值高斯，sigma = 1
%
%   高斯幅度用来模拟组织声阻抗的随机起伏：多数散射子很弱、少数很强，散斑正是
%   这样形成的。随机种子固定，保证 P1、由它派生的囊肿体模（P1 减去囊肿内的
%   散射子）以及参与比较的每一种方法看到的都是同一批散射子——否则方法之间的
%   对比度差异会主要来自体模本身，而不是方法。
%
%   为什么在 x-z 上均匀分布
%   ----------------------
%   均匀 ⇒ 散斑的统计性质与位置无关，这样"同一深度不同横向位置"或"同一横向位置
%   不同深度"的图像差异才只来自成像方法，而不是体模本身的空间起伏。若做成随机
%   聚簇或渐变密度，分辨率与对比度的比较都会被体模结构污染。
%
%   密度：44 x 55 mm 内 20 000 个散射子 = 8.26 个/mm²，折算到一个 −6 dB 分辨
%   单元（0.72 x 0.19 mm）里平均约 1 个（实测 0–4 个）。这是二维体模的已知局限：
%   真实组织每个分辨单元内有远多于 1 个散射子。它对散斑统计与复合 SNR 的影响见
%   报告 5.4 节、5.7 节与附录 B19；用 4 倍密度做对照（`run_dense_snr.m`）得到的
%   结论是"密度不是 SNR 缺额的原因"。

rng(P.randomSeed, 'twister');                 % 固定随机数，与 MATLAB 全局状态无关
% 两个 rand 调用的先后就是随机数的消耗顺序（x 在前、z 在后）。顺序一变，P1 的
% 散射子位置就全变了，所有图像都要重算——改动这里请当作“换体模”。
posX = P.phantomExtentX(1) + diff(P.phantomExtentX) * rand(P.nScatterers, 1);
posZ = P.phantomExtentZ(1) + diff(P.phantomExtentZ) * rand(P.nScatterers, 1);
pos = [posX, zeros(P.nScatterers, 1), posZ];  % y = 0：二维体模
amp = randn(P.nScatterers, 1);                % 零均值高斯幅度（sigma = 1）
end
