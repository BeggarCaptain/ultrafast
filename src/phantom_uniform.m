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
%   密度：44 x 55 mm 内 20 000 个散射子 ≈ 8.3 个/mm²，即每个分辨单元只有
%   几个散射子。这是二维体模的已知局限，报告 5.4 节与
%   docs/实验报告_讨论与附录.md 的 6.5 节讨论了它对散斑统计和复合 SNR 的影响。

rng(P.randomSeed, 'twister');                 % 固定随机数，与 MATLAB 全局状态无关
% 两个 rand 调用的先后就是随机数的消耗顺序（x 在前、z 在后）。顺序一变，P1 的
% 散射子位置就全变了，所有图像都要重算——改动这里请当作“换体模”。
posX = P.phantomExtentX(1) + diff(P.phantomExtentX) * rand(P.nScatterers, 1);
posZ = P.phantomExtentZ(1) + diff(P.phantomExtentZ) * rand(P.nScatterers, 1);
pos = [posX, zeros(P.nScatterers, 1), posZ];  % y = 0：二维体模
amp = randn(P.nScatterers, 1);                % 零均值高斯幅度（sigma = 1）
end
