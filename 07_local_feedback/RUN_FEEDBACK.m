function RUN_FEEDBACK()
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'engine'),fullfile(root,'engine','core'),'-begin');
assert(strcmp(which('CircularPipe'),fullfile(root,'engine','core','CircularPipe.m')), ...
 '存在旧版同名类，请保存工作后clear classes，再运行。');
[u,b]=CONFIG_FEEDBACK();RunFeedbackExperiment(u,b);
end
