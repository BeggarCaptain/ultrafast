function [pos, amp, sel] = phantom_cyst_removed(P, pos0, amp0)
%PHANTOM_CYST_REMOVED  落在九个无回声囊肿内部的散射子。
%
%   [POS, AMP, SEL] = PHANTOM_CYST_REMOVED(P, POS0, AMP0)
%
%   输入
%     P          ：参数结构体（提供囊肿中心与半径）
%     pos0, amp0 ：完整的 P1 散射子集合（由 phantom_uniform 生成）
%
%   输出
%     pos, amp   ：落在囊肿内部的散射子子集
%     sel        ：作用于 P1 散射子的逻辑掩膜（便于统计与核对）
%
%   为什么用子集而不是再做一次全量体模
%   ----------------------------------
%   Field II 对散射子幅度是线性的，所以“P1 上挖九个洞”的体模通道数据恰好等于
%       P2 = P1 - REM
%   其中 REM 是只对“被挖掉的散射子”做一次采集。两次采集（全量 + 被移除部分，
%   各 199 次发射）就同时得到两个体模，省掉第三次全量仿真。相减是在波束形成
%   之后、取包络之前做的（见 stage_analysis），因此严格成立。
%
%   九个囊肿是 P.cystX 与 P.cystZ 的笛卡尔积（3 x 3）；如果只取其中三个
%   （早期版本的一个 bug），所谓“无回声体模”里还剩六个囊肿，对比度测量就废了。

C = cyst_centres(P);                        % 9 x 2 [x z]
sel = false(size(pos0, 1), 1);
for k = 1:size(C, 1)
    d2 = (pos0(:,1) - C(k,1)).^2 + (pos0(:,3) - C(k,2)).^2;   % 平面内距离平方
    sel = sel | (d2 <= P.cystR^2);          % 落在该囊肿内（九个取并集）
end
pos = pos0(sel, :);                         % 需要被移除的散射子
amp = amp0(sel);
end
