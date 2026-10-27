function [pos, amp] = phantom_points(P)
%PHANTOM_POINTS  P3：八个强点散射子，用于分辨率测试。
%
%   [POS, AMP] = PHANTOM_POINTS(P) 返回
%     POS ：8 x 3 矩阵 [x y z]，单位米（y = 0，二维体模）
%     AMP ：8 x 1 幅度列向量（P.pointTargetAmp，取得较强，即使不加噪声 PSF 也远高于
%           数值噪声底）
%
%   八个目标中六个在轴上，覆盖 15–55 mm 的深度；另外两个在 x = ±12 mm、
%   z = 30 mm，用来检查横向均匀性。这正好是 E4 的“分辨率随深度”和“分辨率随
%   横向位置”两条曲线需要的数据。
%   Field II 的计算量随散射子个数增长，这个体模只有 8 个点，几乎不花时间，
%   所以直接按 256 线采集（点目标 PSF 采样更细）。

pos = [P.pointTargetX(:), zeros(numel(P.pointTargetX), 1), P.pointTargetZ(:)];   % [x 0 z]
amp = P.pointTargetAmp * ones(numel(P.pointTargetX), 1);                % 等幅、强幅度
end
