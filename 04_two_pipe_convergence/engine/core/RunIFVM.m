function result = RunIFVM(cfg)

if nargin == 0
    cfg = IFVMConfig();
end
settings=MethodSettings(cfg.method);
mode=settings.method;
validateattributes(cfg.dx, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(cfg.tf, {'numeric'}, {'scalar','real','finite','positive'});
nPool = numel(cfg.L);
nGate = nPool + 1;
if nPool<1, error('IFVM:PoolCount','At least one pool is required.'); end
nCell = round(cfg.L/cfg.dx);
if any(abs(nCell*cfg.dx-cfg.L) > 1e-7) || any(nCell < 2)
    error('IFVM:Grid', 'Each pool length must be an integer multiple of dx, with at least two cells.');
end
if numel(cfg.closingDuration) ~= nGate || numel(cfg.finalOpening) ~= nGate
    error('IFVM:GateSize', 'closingDuration/finalOpening must have nPool+1 entries.');
end

for pool = 1:nPool
    if strcmp(cfg.geometry,'citygate')
        xsect(pool)=CityGate(cfg.sectionD(pool),cfg.sectionH(pool),cfg.slotRatio);
    else
        xsect(pool)=CircularPipe(cfg.sectionD(pool),cfg.slotRatio);
    end
end
A = cell(nPool,1); Q = cell(nPool,1); y = cell(nPool,1);
x = cell(nPool,1); xReal = cell(nPool,1); zReal = cell(nPool,1);
uMax = -inf; aMax = 0;
for pool = 1:nPool
    if isempty(cfg.initialState)
        y{pool}=linspace(cfg.initialDepthUp(pool),cfg.initialDepthDown(pool),nCell(pool));
        Q{pool}=cfg.QDesign(pool)*ones(1,nCell(pool));
    else
        y{pool}=cfg.initialState.y{pool}(:).'; Q{pool}=cfg.initialState.Q{pool}(:).';
        if numel(y{pool})~=nCell(pool) || numel(Q{pool})~=nCell(pool)
            error('IFVM:InitialState','Initial-state grid differs from cfg.dx.');
        end
    end
    A{pool} = xsect(pool).getAofY(y{pool});
    uMax = max(uMax, max(abs(Q{pool})./A{pool}));
    aMax = max(aMax, xsect(pool).a);
    x{pool} = (0.5:1:nCell(pool)-0.5)*cfg.dx;
    xReal{pool} = x{pool} + sum(cfg.L(1:pool-1));
end
for pool = nPool:-1:1

    zReal{pool} = sum(cfg.L(pool+1:end).*cfg.bedSlope(pool+1:end)) ...
        +(cfg.L(pool)-x{pool})*cfg.bedSlope(pool);
end
dt = cfg.dt;
if isempty(dt)
    dt = cfg.Cr*cfg.dx/(uMax+aMax);
end
validateattributes(dt, {'numeric'}, {'scalar','real','finite','positive'});

heads0=externalHeads(0);
[Hu0,Hd0] = GetGateHeads(y, heads0(1), heads0(2), cfg.bedSlope, cfg.dx);
known0=PrescribedFlows(cfg,0,nGate);
G0 = zeros(nGate,1);
targetGateFlow = zeros(nGate,1);
for gate = 1:nGate
    sectionIndex = max(1,min(gate-1,nPool));
    targetQ = cfg.QDesign(sectionIndex);
    targetGateFlow(gate) = targetQ;
    if ~isempty(cfg.gateOpening0)
        G0(gate)=cfg.gateOpening0(gate); continue;
    elseif isfinite(known0(gate))
        G0(gate)=cfg.sectionD(sectionIndex); targetGateFlow(gate)=known0(gate); continue;
    end
    dH = Hu0(gate)-Hd0(gate);
    if strcmp(cfg.initialGateLaw,'classical')
        psi=sign(dH)*sqrt(abs(dH));
    elseif strcmp(cfg.initialGateLaw,'smooth')
        psi=dH/(dH^2+cfg.headEps^2)^(1/4);
    else
        error('IFVM:InitialGateLaw','initialGateLaw must be classical or smooth.');
    end
    requiredArea = abs(targetQ)/(cfg.Cd*sqrt(2*9.8)*abs(psi));
    if targetQ*dH <= 0 || ~isfinite(requiredArea) || requiredArea > xsect(sectionIndex).Af
        error('IFVM:InitialGate', 'Gate %d initial head and design flow are incompatible.',gate);
    end
    G0(gate) = xsect(sectionIndex).getYofA(requiredArea);
end

nStep=ceil(cfg.tf/dt); t=(0:nStep)*dt; t(end)=cfg.tf;
if isfield(cfg,'eventTimes')
    ev=cfg.eventTimes(:).'; ev=ev(ev>0 & ev<cfg.tf); t=[t,ev];
end
t=sort(t); keep=[true,diff(t)>1e-11*max(1,cfg.tf)]; t=t(keep);
if t(end)<cfg.tf, t(end+1)=cfg.tf; end
nStep=numel(t)-1;
nState = nStep+1;
G = zeros(nGate,nState);
for j = 1:nState
    G(:,j) = openingAt(t(j));
end

qBoundary = zeros(nGate,nStep);
qInstant = zeros(nGate,nState);
qOldUsed = zeros(nGate,nStep);
qTarget = zeros(nGate,nStep);
qResidual = zeros(nGate,nStep);
Qofftake = zeros(nPool,nStep);
depthUp = zeros(nPool,nState); depthDown = zeros(nPool,nState);
cellQUp = zeros(nPool,nState); cellQDown = zeros(nPool,nState);
volume = zeros(nPool,nState);
massResidual = zeros(nPool,nStep);
gateHu = zeros(nGate,nState); gateHd = zeros(nGate,nState);
residualNorm = zeros(1,nStep); tolerance = zeros(1,nStep);
iterations = zeros(1,nStep); evaluations = zeros(1,nStep); methodCode = zeros(1,nStep);
if cfg.storeFullField
    yHistory = cell(nPool,1); QHistory = cell(nPool,1);
    for pool=1:nPool
        yHistory{pool}=zeros(nCell(pool),nState);
        QHistory{pool}=zeros(nCell(pool),nState);
    end
else
    yHistory={}; QHistory={};
end
[qInstant(:,1),gateHu(:,1),gateHd(:,1)] = GateFormulaFromState( ...
    y,G(:,1),heads0(1),heads0(2),cfg.bedSlope,cfg.dx, ...
    xsect,cfg.Cd,cfg.headEps,cfg.allowReverse,settings.gateLaw,Q,cfg.gatePhysics);
qInstant(isfinite(known0),1)=known0(isfinite(known0));
recordState(1);
qInitial = qInstant(:,1);
options = cfg.solverOptions;
options.method = mode;
options.gatePhysics=cfg.gatePhysics;
if isfield(cfg,'numerics'),options.numerics=cfg.numerics;end
options.externalHeadOld=heads0;
options.boundaryMomentumDepth=cfg.boundaryMomentumDepth;
if ~cfg.quiet
    fprintf('[%s/%s] dt=%.9g s, tf=%.6g s, steps=%d, dx=%.9g m\n',cfg.caseID,mode,dt,cfg.tf,nStep,cfg.dx);
end
plannedSteps=nStep; completed=0;
failure=struct(); snapshots={}; captured=false(size(cfg.captureTimes));
worstSnapshot=[]; worstIterations=-1;
timer = tic;
for j=1:nStep
    dtStep=t(j+1)-t(j);
    if isempty(cfg.offtakeFunction)
        off=cfg.Qofftake(:);
    else
        off=cfg.offtakeFunction(t(j));
        off=off(:);
    end
    if numel(off)~=nPool || any(~isfinite(off))
        error('IFVM:Offtake','offtakeFunction must return nPool finite flows.');
    end
    options.gateOpeningOld=G(:,j);
    options.externalHeadOld=externalHeads(t(j));
    headsNew=externalHeads(t(j+1));
    options.fixedFlowOld=PrescribedFlows(cfg,t(j),nGate);
    options.fixedFlowNew=PrescribedFlows(cfg,t(j+1),nGate);
    take=any(~captured & cfg.captureTimes<=t(j)+1e-10);
    oldForSnapshot=struct('A',{A},'Q',{Q},'y',{y});
    options.recordConvergence=false;
    if isfield(cfg,'convergenceTarget') && numel(cfg.convergenceTarget)==2
        tt=cfg.convergenceTarget;
        matchTol=1e-8*max([1,abs(tt(:).')]);
        options.recordConvergence=abs(t(j)-tt(1))<=matchTol && abs(t(j+1)-tt(2))<=matchTol;
    end

    try
        [state,stats] = SolveGateCoupling( ...
            qInitial,G(:,j+1),A,Q,y,cfg.nManning,cfg.bedSlope, ...
            cfg.dx,dtStep,xsect,off,cfg.Cd,headsNew(1), ...
            headsNew(2),cfg.headEps,cfg.allowReverse,options);
    catch exception
        failure.step=j; failure.timeStart=t(j); failure.timeEnd=t(j+1);
        failure.message=exception.message; failure.Aold=A; failure.Qold=Q; failure.yold=y;
        break;
    end
    snap=struct('timeStart',t(j),'timeEnd',t(j+1),'old',oldForSnapshot, ...
        'qInitial',qInitial,'openingOld',G(:,j),'opening',G(:,j+1), ...
        'options',options,'headsNew',headsNew,'offtake',off,'stats',stats);
    if options.recordConvergence || ~stats.converged
        result.convergenceSnapshot=snap;
    end
    if take
        snapshots{end+1}=snap; captured(cfg.captureTimes<=t(j)+1e-10)=true;
    end
    if stats.iterations>worstIterations
        worstSnapshot=snap; worstIterations=stats.iterations;
    end
    if ~stats.converged
        failure.step=j; failure.timeStart=t(j); failure.timeEnd=t(j+1);
        failure.message=stats.message; failure.stats=stats;
        failure.Aold=A; failure.Qold=Q; failure.yold=y;
        break;
    end
    A=state.A; Q=state.Q; y=state.y;
    qBoundary(:,j)=state.qGate;
    qInstant(:,j+1)=state.qFormula;
    qOldUsed(:,j)=state.qFormulaOld;
    qTarget(:,j)=state.qTarget;
    qResidual(:,j)=state.residual;
    gateHu(:,j+1)=state.Hu; gateHd(:,j+1)=state.Hd;
    Qofftake(:,j)=off;
    residualNorm(j)=stats.residualNorm;
    tolerance(j)=options.tolAbs+options.tolRel*max(options.flowScale,norm(state.qTarget,inf));
    iterations(j)=stats.iterations; evaluations(j)=stats.evaluations; methodCode(j)=stats.methodCode;
    recordState(j+1);
    massResidual(:,j)=volume(:,j+1)-volume(:,j) ...
        -dtStep*(state.qGate(1:nPool)-state.qGate(2:nGate)-off);

    qInitial=state.qGate;
    completed=j;
    if ~cfg.quiet && cfg.progressEvery>0 && (mod(j,cfg.progressEvery)==0 || j==nStep)
        fprintf('[%s] %d/%d, t=%.6g s, residual=%.3g m^3/s\n', ...
            mode,j,nStep,t(j+1),residualNorm(j));
    end
end
result.computeSeconds=toc(timer);
if completed==plannedSteps, status='completed'; else, status='failed'; end
stateIndex=1:completed+1; stepIndex=1:completed;
t=t(stateIndex); G=G(:,stateIndex);
qInstant=qInstant(:,stateIndex);
depthUp=depthUp(:,stateIndex);
depthDown=depthDown(:,stateIndex);
cellQUp=cellQUp(:,stateIndex);
cellQDown=cellQDown(:,stateIndex);
volume=volume(:,stateIndex);
gateHu=gateHu(:,stateIndex);
gateHd=gateHd(:,stateIndex);
qBoundary=qBoundary(:,stepIndex);
qOldUsed=qOldUsed(:,stepIndex);
qTarget=qTarget(:,stepIndex);
qResidual=qResidual(:,stepIndex);
Qofftake=Qofftake(:,stepIndex);
massResidual=massResidual(:,stepIndex);
residualNorm=residualNorm(:,stepIndex);
tolerance=tolerance(:,stepIndex);
iterations=iterations(:,stepIndex);
evaluations=evaluations(:,stepIndex);
methodCode=methodCode(:,stepIndex);
if cfg.storeFullField
    for pool=1:nPool
        yHistory{pool}=yHistory{pool}(:,stateIndex);
        QHistory{pool}=QHistory{pool}(:,stateIndex);
    end
end
result.status=status; result.failure=failure; result.completedSteps=completed;
result.plannedSteps=plannedSteps; result.method=mode; result.gateLaw=settings.gateLaw;
result.timeCoupling=settings.timeCoupling;
result.requiresClosure=~strcmp(settings.relaxation,'none');

result.version='IFVM_Trajkovic_Study_20261105_v2';
result.snapshots=snapshots; result.worstSnapshot=worstSnapshot;
result.initialMode=cfg.initialMode;
result.config=cfg;
result.dtNominal=dt; result.t=t; result.dt=diff(t);
result.tFlux=0.5*(t(1:end-1)+t(2:end));
result.G0=G0; result.G=G;
result.qBoundary=qBoundary; result.qInstant=qInstant;
result.qOldUsed=qOldUsed; result.qTarget=qTarget; result.qResidual=qResidual;
result.gateHu=gateHu; result.gateHd=gateHd;
result.Qofftake=Qofftake;
result.depthUp=depthUp; result.depthDown=depthDown;
result.cellQUp=cellQUp; result.cellQDown=cellQDown;
result.volume=volume; result.massResidual=massResidual;
result.residualNorm=residualNorm; result.tolerance=tolerance;
result.iterations=iterations; result.evaluations=evaluations; result.methodCode=methodCode;
result.y=yHistory; result.Q=QHistory;
result.finalA=A; result.finalQ=Q; result.finalY=y;
result.x=x; result.xReal=xReal; result.zReal=zReal;
result.initialHeadCompatibilityError=max(abs(qInstant(:,1)-targetGateFlow));
result.initialGateHeads=[Hu0,Hd0];
result.pressureWaveSpeed=zeros(1,nPool);
for pool=1:nPool, result.pressureWaveSpeed(pool)=xsect(pool).a; end
result.pressureTravelTime=cfg.L./result.pressureWaveSpeed;
result.dtOverPressureTravel=dt./result.pressureTravelTime;
if ~cfg.quiet
    fprintf('[%s] %s at t=%.9g s, compute %.3f s, max mass residual %.3g m^3\n', ...
        mode,status,t(end),result.computeSeconds,max([0;abs(massResidual(:))]));
    if strcmp(status,'failed'), fprintf('Reason: %s\n',failure.message); end
end
result.totalEvaluations=sum(evaluations);
if isfield(failure,'stats'), result.totalEvaluations=result.totalEvaluations+failure.stats.evaluations; end

    function h=externalHeads(time)
        h=[cfg.HupExternal;cfg.HdownExternal];
        if ~isempty(cfg.externalHeadFunction), h=cfg.externalHeadFunction(time); h=h(:); end
        if numel(h)~=2 || any(~isfinite(h)), error('IFVM:Head','External heads invalid.'); end
    end

    function opening = openingAt(time)
        if isfield(cfg,'singleOpening') && isfinite(cfg.singleOpening) && ~cfg.singleHold
            opening=TrajOpening(time,G0,cfg.singleOpening,cfg.singleClosureTime, ...
                cfg.reopenTime,cfg.reopeningDuration);
        elseif ~isempty(cfg.openingFunction)
            opening=cfg.openingFunction(time,G0); opening=opening(:);
        elseif isfield(cfg,'operation') && ~isempty(cfg.operation)
            opening=GateSchedule(time,G0,cfg.operation);
        else
            opening=zeros(nGate,1);
            for g=1:nGate
                opening(g)=t2G(G0(g),cfg.finalOpening(g),cfg.closingDuration(g),time);
            end
        end
        if numel(opening)~=nGate || any(~isfinite(opening)) || any(opening<0)
            error('IFVM:Opening','openingFunction must return nPool+1 finite nonnegative openings.');
        end
    end

    function recordState(index)
        for p=1:nPool
            depthUp(p,index)=y{p}(1); depthDown(p,index)=y{p}(end);
            cellQUp(p,index)=Q{p}(1); cellQDown(p,index)=Q{p}(end);
            volume(p,index)=cfg.dx*sum(A{p});
            if cfg.storeFullField
                yHistory{p}(:,index)=y{p}(:);
                QHistory{p}(:,index)=Q{p}(:);
            end
        end
    end
end
