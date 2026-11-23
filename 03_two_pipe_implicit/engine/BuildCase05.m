function c=BuildCase05(u)
validateattributes(u.Q0,{'numeric'},{'scalar','positive','finite'});
validateattributes(u.dx,{'numeric'},{'scalar','positive','finite'});
assert(numel(u.L)==2 && numel(u.D)==2 && numel(u.slope)==2 && numel(u.n)==2,'本案例固定为两根串联圆管。');
assert(all(u.L>0) && all(u.D>0) && all(u.n>0) && all(u.initialGateHeadLoss>0),'管长、管径、糙率和初始闸门损失须为正。');
assert(numel(u.initialGateHeadLoss)==3 && u.endGateStart>0 && all(u.endGateDuration>0),'闸门时间或初始损失设置错误。');
assert(u.middleGateStart>=u.endGateStart && u.middleGateDuration>0,'中闸缓关不能早于初期保持结束。');
assert(u.middleGateFinalFactor>=0 && u.middleGateFinalFactor<=1,'中闸最终开度比例须在0与1之间。');
assert(u.tf>max([u.endGateStart+u.endGateDuration,u.middleGateStart+u.middleGateDuration]),'总时长须覆盖全部动作和最终保持。');
assert(all(abs(u.L/u.dx-round(u.L/u.dx))<1e-9) && all(u.L/u.dx>=2),'管长须为网格长度的整数倍。');
c=BaseConfig('two_strong','screen');c.caseID='case05_two_joint_closure';c.method='E1';
c.geometry='circular';c.gatePhysics='submerged';c.numerics=struct('kernel','rusanov','boundary','riemann');
c.L=u.L;c.sectionD=u.D;c.bedSlope=u.slope;c.nManning=u.n;c.dx=u.dx;
c.slotRatio=u.slotRatio;c.Cr=u.Cr;c.dt=[];c.Cd=u.Cd;c.QDesign=u.Q0*[1,1];
c.Qofftake=[0,0];c.initialState=struct('y',{{}},'Q',{{}});h=u.inletDepth;
for p=1:2
 xs=CircularPipe(u.D(p),u.slotRatio);x=(.5:round(u.L(p)/u.dx)-.5)*u.dx;
 rhs=@(xx,y)(u.slope(p)-xs.getSf(u.n(p),y,u.Q0))/ ...
  (1-(u.Q0/xs.getAofY(y)/xs.getCofY(y))^2);
 opts=odeset('RelTol',1e-9,'AbsTol',1e-11);
 [xx,yy]=ode45(rhs,[0,u.L(p)],h,opts);
 assert(all(yy>u.D(p)),'本工况初始设为有压稳态，请提高inletDepth或减小初始中闸损失。');
 c.initialState.y{p}=interp1(xx,yy,x,'pchip');c.initialState.Q{p}=u.Q0*ones(size(x));
 if p==1,h=yy(end)-u.initialGateHeadLoss(2);end
end
c.initialDepthUp=[c.initialState.y{1}(1),c.initialState.y{2}(1)];
c.initialDepthDown=[c.initialState.y{1}(end),c.initialState.y{2}(end)];
c.HupExternal=c.initialDepthUp(1)-.5*u.dx*u.slope(1)+u.initialGateHeadLoss(1);
c.HdownExternal=c.initialDepthDown(2)+.5*u.dx*u.slope(2)-u.initialGateHeadLoss(3);
c.tf=u.tf;c.operation=[];c.captureTimes=[];c.storeFullField=false;
c.initialMode='discrete_steady_at_given_flow';c.initialGateLaw='classical';c.gateOpening0=[];
c.prescribedFlowConstant=nan(3,1);c.boundaryMomentumDepth=[nan,nan];
c.progressEvery=u.progressEvery;c.quiet=false;
c.closingDuration=[inf,inf,inf];c.finalOpening=[0,0,0];
c.openingFunction=@(t,G0)OpeningCase05(t,G0,u);
c.eventTimes=unique([u.endGateStart,u.endGateStart+u.endGateDuration,u.middleGateStart,u.middleGateStart+u.middleGateDuration]);
end
