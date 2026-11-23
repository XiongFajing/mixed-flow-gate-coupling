function folder=RunSmoothingBatch(u,b)
validateattributes(b.multipliers,{'numeric'},{'vector','positive','finite','increasing'});
validateattributes(b.baseHeadEps,{'numeric'},{'scalar','positive','finite'});
assert(any(strcmp(b.method,{'S1','S2'})),'method仅支持S1或S2。');
assert(numel(b.errorWindow)==2 && b.errorWindow(1)>=0 && ...
 b.errorWindow(2)<=u.tf && diff(b.errorWindow)>0,'误差时段错误。');
if b.referenceEnabled
 assert(numel(b.referenceHeadEps)==2 && all(b.referenceHeadEps>0) && ...
 b.referenceHeadEps(2)<b.referenceHeadEps(1),'参考eps须为两个递减正数。');
 assert(min(b.baseHeadEps*b.multipliers)>= ...
 b.referenceMinEpsRatio*b.referenceHeadEps(2),'参考eps不够小，请调整参考参数。');
end
if ~isfolder(b.outputRoot),mkdir(b.outputRoot);end
key=struct('u',u,'b',b); key.b=rmfield(key.b,{'resume','rerunFailed','forceSteady','outputRoot'});
folder='';
if b.resume
 d=dir(fullfile(b.outputRoot,'batch_*'));d=d([d.isdir]);
 [~,order]=sort([d.datenum],'descend');
 for j=order
  candidate=fullfile(d(j).folder,d(j).name); f=fullfile(candidate,'batch_config.mat');
  if isfile(f)
   s=load(f,'key');
   if isfield(s,'key') && isequaln(s.key,key),folder=candidate;break;end
  end
 end
end
if isempty(folder)
 folder=tempname(b.outputRoot);[~,name]=fileparts(folder);
 folder=fullfile(b.outputRoot,['batch_',datestr(now,'yyyymmdd_HHMMSS'),'_',name]);mkdir(folder);
end
save(fullfile(folder,'batch_config.mat'),'u','b','key');
engine=fileparts(mfilename('fullpath')); prepTimer=tic;
cache=fullfile(engine,'initial_cache.mat');
prep=u;prep.Cr=u.steadyCr;

steadyKey=rmfield(prep,{'solver','headEps','makeFigures','figureVisible','progressEvery'});
hit=false;
if ~b.forceSteady && isfile(cache)
 s=load(cache);
 if isequaln(s.steadyKey,steadyKey),c=s.c;info=s.info;hit=true;end
end
if ~hit
 c=BuildCase05(prep);[c,info]=PrepareCase05Steady(c,prep);
 save(cache,'steadyKey','c','info');
end
prepareSeconds=toc(prepTimer);
c.Cr=u.Cr;c.dt=[];c.quiet=true;c.captureTimes=[];c.storeFullField=false;
fn=fieldnames(u.solver);
for j=1:numel(fn),c.solverOptions.(fn{j})=u.solver.(fn{j});end
Href=abs(info.initialGateHeadDifference(2));Qref=abs(u.Q0);
assert(Href>0 && Qref>0,'参考尺度须为正。');

vmax=0;amax=0;
for p=1:2
 xs=CircularPipe(c.sectionD(p),c.slotRatio);
 vmax=max(vmax,max(abs(c.initialState.Q{p})./xs.getAofY(c.initialState.y{p})));
 amax=max(amax,xs.a);
end
c.dt=c.Cr*c.dx/(vmax+amax);
save(fullfile(folder,'initial_state.mat'),'c','info','Href','Qref','prepareSeconds','hit');
fprintf('批次目录：%s\nHref=%.10g m；Qref=%.10g m3/s；dt=%.10g s\n',folder,Href,Qref,c.dt);

for j=1:numel(b.multipliers)
 epsH=b.baseHeadEps*b.multipliers(j);
 cc=c;cc.method=b.method;cc.headEps=epsH;
 group=fullfile(folder,sprintf('group_%02d',j));
 RunGroup(cc,group,b,epsH/Href,b.multipliers(j),'批量');
 SummarizeSmoothing(folder);
end
if b.referenceEnabled
 for j=1:2
  cc=c;cc.method='S2';cc.headEps=b.referenceHeadEps(j);
  cc.solverOptions.tolAbs=b.referenceTolAbs;
  cc.solverOptions.tolRel=b.referenceTolRel;
  cc.solverOptions.picardMaxIter=b.referenceMaxIter;
  RunGroup(cc,fullfile(folder,sprintf('reference_%d',j)),b,cc.headEps/Href,NaN,'参考');
 end
end
SummarizeSmoothing(folder);
fprintf('批量结束；请查看 metrics.csv 和 reference_check.csv。\n');
end

function RunGroup(c,folder,b,alpha,multiplier,kind)
f=fullfile(folder,'results.mat');
if b.resume && isfile(f)
 old=load(f,'meta');
 if ~(b.rerunFailed && ~strcmp(old.meta.status,'completed'))
  fprintf('复用已有组：%s\n',folder);return;
 end
end
if ~isfolder(folder),mkdir(folder);end
fprintf('%s：%s，epsilon=%.8g m，alpha=%.8g\n',kind,c.method,c.headEps,alpha);
meta=struct('status','exception','message','','alpha',alpha,'multiplier',multiplier, ...
 'headEps',c.headEps,'method',c.method,'formalSeconds',NaN,'exportSeconds',NaN);
r=[];timer=tic;
try
 r=RunIFVM(c);meta.formalSeconds=toc(timer);meta.status=r.status;
 if strcmp(r.status,'failed'),meta.message=r.failure.message;end
catch ex
 meta.formalSeconds=toc(timer);meta.message=getReport(ex,'extended','hyperlinks','off');
end
timer=tic;

save(f,'r','meta','c','-v7.3');
meta.exportSeconds=toc(timer);save(f,'meta','-append');
fprintf('状态%s，正式耗时%.4g s\n',meta.status,meta.formalSeconds);
end
