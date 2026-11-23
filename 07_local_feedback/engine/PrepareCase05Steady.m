function [c,info]=PrepareCase05Steady(c,u)

assert(u.steadyChunk>0 && u.steadyMaxTime>=2*u.steadyChunk,'稳态准备至少需要两个检查时间段。');
info.elapsed=zeros(1,2);info.headRate=zeros(1,2);info.flowRate=zeros(1,2);
for p=1:2
 w=c;w.L=c.L(p);w.sectionD=c.sectionD(p);w.bedSlope=c.bedSlope(p);w.nManning=c.nManning(p);
 w.QDesign=u.Q0;w.Qofftake=0;w.prescribedFlowConstant=[u.Q0;u.Q0];
 w.closingDuration=[inf,inf];w.finalOpening=[0,0];w.gateOpening0=[c.sectionD(p);c.sectionD(p)];
 w.initialDepthUp=c.initialDepthUp(p);w.initialDepthDown=c.initialDepthDown(p);
 w.initialState.y={c.initialState.y{p}};w.initialState.Q={c.initialState.Q{p}};
 w.openingFunction=[];w.operation.times=[0,u.steadyChunk];w.operation.factors=ones(2,2);
 w.tf=u.steadyChunk;w.eventTimes=[];w.quiet=true;w.captureTimes=[];w.storeFullField=false;
 w.HupExternal=w.initialState.y{1}(1)+.04;w.HdownExternal=w.initialState.y{1}(end)-.04;
 converged=false;previousPassed=false;total=0;
 while total+u.steadyChunk<=u.steadyMaxTime+1e-10
  y0=w.initialState.y{1};q0=w.initialState.Q{1};r=RunIFVM(w);
  if ~strcmp(r.status,'completed'),error('管道%d稳态准备失败：%s',p,r.failure.message);end

  dy=max(abs(r.finalY{1}-y0))/u.steadyChunk;
  dq=max(abs(r.finalQ{1}-q0))/u.steadyChunk;
  total=total+u.steadyChunk;w.initialState.y=r.finalY;w.initialState.Q=r.finalQ;
  passed=dy<=u.steadyHeadRateTol && dq<=u.steadyFlowRateTol;
  if passed && previousPassed,converged=true;break;end
  previousPassed=passed;
 end
 if ~converged
  error('管道%d在%.3g s准备时间内未达到稳态：水头变化速率%.3g m/s，流量变化速率%.3g m³/s²。可增大steadyMaxTime。',p,total,dy,dq);
 end
 c.initialState.y{p}=w.initialState.y{1};c.initialState.Q{p}=w.initialState.Q{1};
 info.elapsed(p)=total;info.headRate(p)=dy;info.flowRate(p)=dq;
 fprintf('管道%d稳态准备完成，准备时长%.3g s。\n',p,total);
end
c.initialDepthUp=[c.initialState.y{1}(1),c.initialState.y{2}(1)];
c.initialDepthDown=[c.initialState.y{1}(end),c.initialState.y{2}(end)];

c.HupExternal=c.initialDepthUp(1)-.5*c.dx*c.bedSlope(1)+u.initialGateHeadLoss(1);
c.HdownExternal=c.initialDepthDown(2)+.5*c.dx*c.bedSlope(2)-u.initialGateHeadLoss(3);
[hu,hd]=GetGateHeads(c.initialState.y,c.HupExternal,c.HdownExternal,c.bedSlope,c.dx);
c.gateOpening0=zeros(3,1);
for g=1:3
 p=max(1,min(g-1,2));xs=CircularPipe(c.sectionD(p),c.slotRatio);
 dh=hu(g)-hd(g);assert(dh>0,'稳态中闸头差非正，请调整初始水位与中闸损失。');
 area=u.Q0/(c.Cd*sqrt(2*9.8*dh));
 assert(area<=xs.Af,'闸门所需开度超过管径，请增大初始头差。');
 c.gateOpening0(g)=xs.getYofA(area);
end
info.initialGateHeadDifference=hu-hd;
end
