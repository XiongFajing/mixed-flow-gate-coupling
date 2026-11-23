function folder=RunFeedbackExperiment(u,b)
assert(u.Cr==.5,'用于捕获共用状态的基准轨迹Cr须为0.5。');
validateattributes(b.timeStepFactors,{'numeric'},{'vector','positive','finite','<=',1});
validateattributes(b.targetTime,{'numeric'},{'scalar','positive','finite','<',u.tf});
validateattributes(b.perturbationRatios,{'numeric'},{'vector','positive','finite'});

if ~isfolder(b.outputRoot),mkdir(b.outputRoot);end
folder=b.outputRoot;
if ~isfolder(folder),mkdir(folder);end
key=struct('u',u,'targetTime',b.targetTime,'version',1);
cache=fullfile(b.outputRoot,'state_cache.mat');hit=false;prepTimer=tic;
if ~b.rebuildTrajectory && isfile(cache)
 s=load(cache);
 if isfield(s,'key') && isequaln(s.key,key)
  c=s.c;snap=s.snap;info=s.info;trace=s.trace;hit=true;
 end
end
if ~hit
 c=BuildCase05(u);
 [c,info]=PrepareCase05Steady(c,u);
 c.method='E1';c.Cr=.5;c.dt=[];c.quiet=false;c.storeFullField=false;

 vmax=0;amax=0;
 for p=1:2
  xs=CircularPipe(c.sectionD(p),c.slotRatio);
  vmax=max(vmax,max(abs(c.initialState.Q{p})./xs.getAofY(c.initialState.y{p})));
  amax=max(amax,xs.a);
 end
 c.dt=c.Cr*c.dx/(vmax+amax);
 c.captureTimes=b.targetTime;

 c.tf=min(u.tf,b.targetTime+3*c.dt);
 r=RunIFVM(c);
 assert(~isempty(r.snapshots),'未捕获目标状态：%s',FailureMessage(r));
 snap=r.snapshots{1};
 assert(snap.timeStart>=b.targetTime-1e-9 && snap.timeEnd>snap.timeStart,'状态捕获错误。');
 trace=struct('t',r.t,'G',r.G,'qBoundary',r.qBoundary, ...
  'gateHu',r.gateHu,'gateHd',r.gateHd,'status',r.status,'dtNominal',r.dtNominal);
 save(cache,'key','c','snap','info','trace','-v7.3');
end
prepareSeconds=toc(prepTimer);
dt0=snap.timeEnd-snap.timeStart;
dt=dt0;
assert(abs(dt-c.dt)<1e-9*max(1,c.dt), ...
 '目标步恰好位于事件/末端短步，请改变targetTime。');
qBase=snap.stats.lastQ(:);opening=snap.opening(:);
assert(opening(2)>0,'所选时刻中闸已关闭，请改变targetTime。');
for p=1:2,xsect(p)=CircularPipe(c.sectionD(p),c.slotRatio);end
[qCheck,huOld,hdOld]=GateFormulaFromState(snap.old.y,opening, ...
 snap.options.externalHeadOld(1),snap.options.externalHeadOld(2), ...
 c.bedSlope,c.dx,xsect,c.Cd,c.headEps,c.allowReverse,'classical',snap.old.Q,c.gatePhysics);
assert(max(abs(qCheck-qBase))<1e-12,'快照流量不是旧头差/步末开度显式流量。');
fprintf('保存状态t=%.12g s；dt=%.12g s；Cr=%.3g\n',snap.timeStart,dt,c.Cr);
fprintf('中闸旧头差=%.9g m；基准流量=%.9g m3/s；开度=%.9g m\n', ...
 huOld(2)-hdOld(2),qBase(2),opening(2));
allSummary=[];allRaw=[];allChecks=[];allTrials=cell(numel(b.timeStepFactors),1);
for d=1:numel(b.timeStepFactors)
 factor=b.timeStepFactors(d);dt=dt0*factor;

 fprintf('核验dt/dt0=%.4g，dt=%.12g s，等效Cr=%.4g\n',factor,dt,c.Cr*factor);
context=struct('gateLaw','classical','gatePhysics',c.gatePhysics, ...
 'qOld',qBase,'boundaryMomentumDepth',snap.options.boundaryMomentumDepth, ...
 'numerics',snap.options.numerics);
context.trialCache=BuildIFVMTrialCache(snap.old.y,snap.old.Q,c.nManning, ...
 c.bedSlope,c.dx,dt,xsect,c.numerics);
nAmplitude=numel(b.perturbationRatios);raw=nan(2*nAmplitude,15);summary=nan(nAmplitude,20);valid=false(2*nAmplitude,1);messages=strings(2*nAmplitude,1);trials=cell(nAmplitude,2);
for a=1:nAmplitude
 dq=b.perturbationRatios(a)*abs(u.Q0);
 for side=1:2
  sg=2*side-3;q=qBase;q(2)=q(2)+sg*dq;timer=tic;
  trial=EvaluateGateCoupling(q,opening,snap.old.A,snap.old.Q,snap.old.y, ...
   c.nManning,c.bedSlope,c.dx,dt,xsect,snap.offtake,c.Cd, ...
   snap.headsNew(1),snap.headsNew(2),c.headEps,c.allowReverse,context);
  seconds=toc(timer);trials{a,side}=trial;k=(a-1)*2+side;
  valid(k)=trial.valid;messages(k)=string(trial.message);mass=[NaN NaN];
  if trial.valid
   for p=1:2
    mass(p)=c.dx*sum(trial.A{p}-snap.old.A{p}) ...
     -dt*(q(p)-q(p+1)-snap.offtake(p));
   end
  end
  raw(k,:)=[b.perturbationRatios(a),dq,sg,snap.timeStart,snap.timeStart+dt, ...
   q.',trial.Hu(2),trial.Hd(2),trial.Hu(2)-trial.Hd(2),trial.qFormula(2),mass,seconds];
 end
 minus=trials{a,1};plus=trials{a,2};
 khu=NaN;khd=NaN;kh=NaN;kq=NaN;slope=NaN;cross=NaN;
 if minus.valid && plus.valid
  khu=(plus.Hu(2)-minus.Hu(2))/(2*dq);
  khd=(plus.Hd(2)-minus.Hd(2))/(2*dq);kh=khu-khd;
  kq=(plus.qFormula(2)-minus.qFormula(2))/(2*dq);
  hm=minus.Hu(2)-minus.Hd(2);hp=plus.Hu(2)-plus.Hd(2);
  if hp~=hm,slope=(plus.qFormula(2)-minus.qFormula(2))/(hp-hm);end
  cross=(min(hm,hp)<=0 && max(hm,hp)>=0);
 else,hm=NaN;hp=NaN;
 end
 summary(a,:)=[b.perturbationRatios(a),dq,snap.timeStart,dt,c.Cr*factor,qBase(2), ...
  huOld(2),hdOld(2),huOld(2)-hdOld(2),khu,khd,kh,kq,slope,hm,hp, ...
  minus.qFormula(2),plus.qFormula(2),cross,minus.valid && plus.valid];
end
assert(all(raw(:,6)==qBase(1)) && all(raw(:,8)==qBase(3)), ...
 '外侧闸门流量未保持一致。');
allSummary=[allSummary;[repmat(factor,size(summary,1),1),summary]];
allRaw=[allRaw;[repmat(factor,size(raw,1),1),raw]];
end
paper=allSummary(abs(allSummary(:,2)-.01)<1e-12,:);
T=array2table(paper(:,[1 4 5 6 2 3 11 12 13 14]),'VariableNames', ...
 {'StepFactor','SnapshotTime_s','Dt_s','EquivalentCr','PerturbationRatio','DeltaQ_m3s','UpstreamResponse_sm2','DownstreamResponse_sm2','HeadDifferenceResponse_sm2','DischargeGain'});
writetable(T,fullfile(folder,'coefficients.csv'));
rawPaper=allRaw(abs(allRaw(:,2)-.01)<1e-12,:);
T=array2table(rawPaper,'VariableNames',{'StepFactor','PerturbationRatio','DeltaQ_m3s','Direction','StartTime_s','EndTime_s','Gate1_Flow_m3s','Gate2_Flow_m3s','Gate3_Flow_m3s','Gate2_UpstreamHead_m','Gate2_DownstreamHead_m','Gate2_HeadDifference_m','Gate2_TargetFlow_m3s','Pool1_MassResidual_m3','Pool2_MassResidual_m3','TrialTime_s'});
writetable(T,fullfile(folder,'trials.csv'));
save(fullfile(folder,'results.mat'),'paper','rawPaper','dt0');
fprintf('输出目录：%s\n',folder);
end

function m=FailureMessage(r)
if isfield(r.failure,'message'),m=r.failure.message;else,m=r.status;end
end
