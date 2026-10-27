function f = stamp_path(P)
%STAMP_PATH  参数指纹文件的路径。
f = fullfile(P.dataDir, 'stamp_acq.mat');
end
