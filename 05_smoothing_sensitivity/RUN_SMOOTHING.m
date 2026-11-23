function RUN_SMOOTHING()
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'engine'),fullfile(root,'engine','core'),'-begin');
[u,b]=CONFIG_SMOOTHING();methods={'S1','S2'};labels={'I2','I3'};
for k=1:2
 b.method=methods{k};b.outputRoot=fullfile(root,'results',labels{k});
 RunSmoothingBatch(u,b);
end
end
