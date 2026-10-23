function C = cyst_centres(P)
%CYST_CENTRES  P.cystX 与 P.cystZ 的笛卡尔积：九个囊肿中心。
%
%   C = CYST_CENTRES(P) 返回 9 x 2 矩阵 [x z]（单位：米）。
%
%   体模定义里把“三列囊肿”（P.cystX）与“三行囊肿”（P.cystZ）分开存放，
%   实际体模在每一个组合处都有一个囊肿，即 3 x 3 = 9 个。散射子剔除
%   （phantom_cyst_removed）与对比度指标（metrics_contrast）必须使用同一组
%   九个中心，因此在这里统一计算一次，避免两边各算一遍、算错一边。
%
%   历史教训：早期版本把 cystX(k) 与 cystZ(k) 一一配对，只生成了 3 个囊肿，
%   对比度体模因此是错的（而且不报错，很难发现）。

[XX, ZZ] = meshgrid(P.cystX, P.cystZ);   % 所有组合，X 沿列变化
C = [XX(:), ZZ(:)];                       % 9 x 2，按行排列
end
