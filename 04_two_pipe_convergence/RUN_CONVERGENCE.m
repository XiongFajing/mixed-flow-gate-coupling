function RUN_CONVERGENCE()
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'engine'),fullfile(root,'engine','core'),'-begin');
u=CONFIG_CONVERGENCE();prep=u;prep.Cr=u.steadyCr;
base=BuildCase05(prep);[base,~]=PrepareCase05Steady(base,prep);
methods={'R1','S1','S2'};labels={'I1','I2','I3'};configs=cell(1,3);
for k=1:3
 c=base;c.method=methods{k};c.Cr=u.Cr;c.dt=[];c.headEps=u.headEps;
 fields=fieldnames(u.solver);
 for j=1:numel(fields),c.solverOptions.(fields{j})=u.solver.(fields{j});end
 c.quiet=true;c.progressEvery=0;c.captureTimes=[];c.storeFullField=false;
 c.convergenceTarget=[];configs{k}=c;
end
folder=fullfile(root,'results');if ~isfolder(folder),mkdir(folder);end
r=RunIFVM(configs{1});
assert(strcmp(r.status,'failed') && isfield(r,'convergenceSnapshot'),'I1 did not produce an iteration failure.');
snap=r.convergenceSnapshot;target=[snap.timeStart,snap.timeEnd];
opts=snap.options;opts.recordConvergence=true;c=configs{1};
for p=1:numel(c.L),xs(p)=CircularPipe(c.sectionD(p),c.slotRatio);end
[~,st]=SolveGateCoupling(snap.qInitial,snap.opening,snap.old.A,snap.old.Q,snap.old.y, ...
 c.nManning,c.bedSlope,c.dx,diff(target),xs,snap.offtake,c.Cd,snap.headsNew(1),snap.headsNew(2),c.headEps,c.allowReverse,opts);
assert(~st.converged,'I1 failure was not reproduced.');
allStats=cell(1,3);allStats{1}=st;
for k=2:3
 c=configs{k};c.tf=target(2);c.convergenceTarget=target;r=RunIFVM(c);
 assert(isfield(r,'convergenceSnapshot'),'Target interval was not reached.');
 snap=r.convergenceSnapshot;
 assert(max(abs([snap.timeStart,snap.timeEnd]-target))<=1e-8,'The physical time intervals differ.');
 allStats{k}=snap.stats;
end
summary=[];
for k=1:3
 st=allStats{k};a=st.convergenceTrace;
 data=[a.iteration,a.residual,a.tolerance,a.normalized,a.dH(:,2),a.Hu(:,2),a.Hd(:,2),a.qTrial(:,2),a.qTarget(:,2),a.gateResidual(:,2),a.omega];
 names={'Evaluation','Residual_m3s','Tolerance_m3s','NormalizedResidual','Gate2_HeadDifference_m','Gate2_UpstreamHead_m','Gate2_DownstreamHead_m','Gate2_TrialFlow_m3s','Gate2_TargetFlow_m3s','Gate2_SignedResidual_m3s','RelaxationFactor'};
 writetable(array2table(data,'VariableNames',names),fullfile(folder,[labels{k},'.csv']));
 row=table(string(labels{k}),target(1),target(2),numel(a.iteration),logical(st.converged),a.residual(end),a.tolerance(end), ...
 'VariableNames',{'Method','StartTime_s','EndTime_s','Evaluations','Converged','Residual_m3s','Tolerance_m3s'});
 if k==1,summary=row;else,summary=[summary;row];end
end
writetable(summary,fullfile(folder,'summary.csv'));
save(fullfile(folder,'results.mat'),'target','allStats');
disp(summary);
end
