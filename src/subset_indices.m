function idx = subset_indices(P)
%SUBSET_INDICES  参与比较的两种方法各自的角度下标。
%
%   IDX = SUBSET_INDICES(P) 返回 1 x numel(P.compoundCounts) 的元胞数组：
%     idx{1} → 单平面波：只取 0° 那一次发射
%     idx{2} → C11：全部 11 个角度
%
%   本实验只比较“单平面波”与“11 角度相干复合”，不再保留 5/9/17/… 等其他
%   档位，所以这里只有两条分支。角度集合由 config/params.m 给出：
%   11 个角度，-10° ~ +10°，间隔 2°。
%
%   为什么只比较"单平面波"与"11 角度复合"两端
%   ------------------------------------------
%   这两端正好是"帧率最高、图像最差"与"帧率降到 1/11、图像最好"的对照，中间
%   档位（5、9 个角度）只是在这条曲线上插值，不提供新的结论。论文里也有多档
%   对比，本实验按需要只保留两端（用户要求不要"随 n 变化"的曲线）。
%
%   为什么单平面波固定取 0°
%   ----------------------
%   0° 是垂直入射，阵元指向性与接收孔径的响应都对称，代表单平面波能达到的最好
%   情况；若取一个倾斜角度，偏转损失会让单平面波显得更差，比较就不公平了。
%
%   为什么用"角度序号"而不是角度值寻址
%   --------------------------------
%   采集文件按序号命名（a001…a011），序号是唯一稳定的键；角度值是浮点数，用它
%   做索引既慢又容易因为舍入而对不上。序号与角度的对应由 params 保证：
%   1 → -10°、6 → 0°、11 → +10°。

I = cell(1, numel(P.compoundCounts));
for k = 1:numel(P.compoundCounts)
    switch P.compoundCounts(k)
        case 1
            I{k} = P.zeroAngleIndex;                  % 单平面波：0° 那一次
        case P.nAngles
            I{k} = 1:P.nAngles;              % 全部角度：C11
        otherwise
            error('subset_indices:unsupported', ...
                '未定义的档位 n=%d（本实验只有 1 与 %d）', P.compoundCounts(k), P.nAngles);
    end
end
idx = I;
end
