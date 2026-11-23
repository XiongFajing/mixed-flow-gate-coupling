function SaveResults(r,folder,formalSeconds)
if ~isfolder(folder),mkdir(folder);end
t=r.t(:);data=t;names={'Time_s'};
for p=1:size(r.depthUp,1)
 data=[data,r.depthUp(p,:).',r.depthDown(p,:).',r.cellQUp(p,:).',r.cellQDown(p,:).'];
 names=[names,{sprintf('Pool%d_UpstreamHead_m',p),sprintf('Pool%d_DownstreamHead_m',p),sprintf('Pool%d_UpstreamCellFlow_m3s',p),sprintf('Pool%d_DownstreamCellFlow_m3s',p)}];
end
for g=1:size(r.G,1)
 data=[data,r.G(g,:).',r.gateHu(g,:).',r.gateHd(g,:).'];
 names=[names,{sprintf('Gate%d_Opening_m',g),sprintf('Gate%d_UpstreamHead_m',g),sprintf('Gate%d_DownstreamHead_m',g)}];
end
if isfield(r,'plotP7') && ~isempty(r.plotP7)
 data=[data,r.plotP7(:)];names=[names,{'P7_Head_m'}];
end
writetable(array2table(data,'VariableNames',names),fullfile(folder,'states.csv'));
names={'StartTime_s','EndTime_s'};
for g=1:size(r.qBoundary,1),names{end+1}=sprintf('Gate%d_Flow_m3s',g);end
data=[r.t(1:end-1).',r.t(2:end).',r.qBoundary.'];
writetable(array2table(data,'VariableNames',names),fullfile(folder,'gate_flows.csv'));
message="";
if isfield(r.failure,'message'),message=string(r.failure.message);end
av=NaN;mx=NaN;
if ~isempty(r.evaluations),av=mean(r.evaluations);mx=max(r.evaluations);end
T=table(string(r.method),string(r.status),r.t(end),r.dtNominal,r.completedSteps,av,mx,formalSeconds,message, ...
 'VariableNames',{'Method','Status','EndTime_s','NominalDt_s','CompletedSteps','MeanEvaluations','MaxEvaluations','FormalTime_s','Message'});
writetable(T,fullfile(folder,'summary.csv'));
drop=intersect(fieldnames(r),{'snapshots','worstSnapshot','y','Q'});
r=rmfield(r,drop);
save(fullfile(folder,'results.mat'),'r','formalSeconds','-v7.3');
end
