function f = stamp_path(P)
%STAMP_PATH  参数指纹文件的路径（data/stamp_acq.mat）。
%
%   单独做成函数，是为了让写指纹的 stamp_guard、冷启动时删指纹的
%   pipeline_reset、以及任何检查都引用同一个路径。否则一旦有人改了文件名，
%   就会出现"写在新名字、删在旧名字"这种两头都对不上的情况——指纹文件还在，
%   数据却已经删光，下次运行会误判为"数据来历不明"。
f = fullfile(P.dataDir, 'stamp_acq.mat');
end
