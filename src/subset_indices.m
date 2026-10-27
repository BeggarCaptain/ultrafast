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
%   注意：索引是“角度序号”，1 对应 -10°、6 对应 0°、11 对应 +10°。
%   采集文件里 a001…a011 就是用这些序号命名的，两边必须一致。

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
