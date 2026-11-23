function folder=RUN_SINGLE()

here=fileparts(mfilename('fullpath'));
addpath(here,fullfile(here,'core'),'-begin');
assert(strcmpi(which('CircularPipe'),fullfile(here,'core','CircularPipe.m')), ...
 '加载了其他包的CircularPipe。请移除旧包路径并clear classes后重试。');
c=CONFIG_SINGLE();
assert(all(isfinite([c.Cd,c.dx,c.Cr,c.tf,c.initialDepth])) && ...
 all([c.Cd,c.dx,c.Cr,c.tf,c.initialDepth]>0),'参数应为有限正数。');
c.caseID='single_A2';c.method='E1';c.geometry='circular';c.dt=[];
c.initialDepthUp=c.initialDepth;c.initialDepthDown=c.initialDepth;
c.initialState=[];c.initialMode='原包初始均匀水深';
c.gateOpening0=[0.014;0.1];c.Qofftake=0;
c.HupExternal=0.1;c.HdownExternal=0;
c.prescribedFlowConstant=[c.QDesign;nan];
c.boundaryMomentumDepth=[c.initialDepth,nan];
c.gatePhysics='free_outlet';c.headEps=0;c.allowReverse=true;
c.closingDuration=[inf,inf];c.finalOpening=[0,0];
c.offtakeFunction=[];c.externalHeadFunction=[];
c.captureTimes=[];c.storeFullField=true;c.quiet=false;c.progressEvery=2000;
c.eventTimes=unique([0,c.singleClosureTime,c.reopenTime, ...
 c.reopenTime+c.reopeningDuration]);
c.numerics=struct('kernel','rusanov','boundary','riemann');

c.solverOptions=struct('tolAbs',1e-6,'tolRel',1e-4, ...
 'flowScale',abs(c.QDesign));
fprintf('A2单管显式方案：下游Cd=%.8f，dx=%.6g m，窄缝比例=%.6g。\n', ...
 c.Cd,c.dx,c.slotRatio);
timer=tic;r=RunIFVM(c);formalSeconds=toc(timer);
r.plotP7=interp1(r.x{1},r.y{1},9.4,'linear');
r.y={};r.Q={};r.snapshots={};r.worstSnapshot=[];
folder=fullfile(here,'results');
mkdir(folder);
SaveResults(r,folder,formalSeconds);
ref=readmatrix(fullfile(here,'reference','A2_P7.csv'));
ok=all(isfinite(ref(:,1:2)),2) & ref(:,1)>=r.t(1) & ref(:,1)<=r.t(end);
obs=ref(ok,1:2);pred=interp1(r.t,r.plotP7,obs(:,1),'linear');pred=pred(:);
err=pred-obs(:,2);[emax,idx]=max(abs(err));
T=table(numel(err),min(obs(:,1)),max(obs(:,1)),sqrt(mean(err.^2)),mean(abs(err)),emax,obs(idx,1),c.Cd, ...
 'VariableNames',{'N','StartTime_s','EndTime_s','RMSE_m','MAE_m','Emax_m','EmaxTime_s','Cd'});
writetable(T,fullfile(folder,'metrics.csv'));
writetable(array2table([obs,pred,err],'VariableNames',{'Time_s','MeasuredHead_m','SimulatedHead_m','Error_m'}),fullfile(folder,'comparison.csv'));
fprintf('结果目录：%s\n',folder);
end
