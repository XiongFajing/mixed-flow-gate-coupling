function cfg=BaseConfig(caseID,profile)

if nargin<1, caseID='two_strong'; end
if nargin<2, profile='screen'; end
cfg.caseID=caseID; cfg.profile=profile; cfg.method='S2';
cfg.geometry='circular'; cfg.L=[10,10]; cfg.sectionD=[0.1,0.1];
cfg.nManning=[0.008,0.008]; cfg.slotRatio=0.01;
cfg.dx=0.5; cfg.Cr=0.5; cfg.dt=[];
cfg.initialDepthUp=[0.11,0.095]; cfg.initialDepthDown=cfg.initialDepthUp;
cfg.QDesign=[0.0013,0.0013]; cfg.initialState=[]; cfg.gateOpening0=[];
cfg.initialMode='uniform_equilibrium'; cfg.initialGateLaw='classical';
cfg.HupExternal=0.15; cfg.HdownExternal=0.055;
cfg.externalHeadFunction=[]; cfg.prescribedFlowFunction=[];
cfg.prescribedFlowConstant=nan(3,1); cfg.prescribedFlowTable=[];
cfg.Cd=0.83; cfg.headEps=2e-4; cfg.allowReverse=true;
cfg.closingDuration=[inf,inf,inf]; cfg.finalOpening=[0,0,0];
cfg.Qofftake=[0,0]; cfg.openingFunction=[]; cfg.offtakeFunction=[];
cfg.operation.times=[0,2,7,12,17,30];
cfg.operation.factors=[ones(2,6);1,1,0.5,0.5,1,1];
cfg.eventTimes=cfg.operation.times;
cfg.tf=30; cfg.storeFullField=true; cfg.progressEvery=200;
cfg.captureTimes=[0,2,3,5,7,12,17]; cfg.quiet=false;
cfg.boundaryMomentumDepth=[nan,nan];
cfg.measurements=[]; cfg.cacheTag='';
cfg.solverOptions.flowScale=max(abs(cfg.QDesign));
cfg.solverOptions.tolAbs=1e-6; cfg.solverOptions.tolRel=1e-4;
cfg.solverOptions.closedOpeningTol=1e-12;
cfg.solverOptions.picardMaxIter=500;
cfg.solverOptions.picardFixedOmega=0.02;
cfg.solverOptions.picardOmegaInitial=0.02;
cfg.solverOptions.picardOmegaMin=1e-4;
cfg.solverOptions.picardOmegaMax=0.5;
cfg.solverOptions.picardMaxFlowStep=0.5*max(abs(cfg.QDesign));
switch caseID
    case 'two_normal'
    case 'two_strong'
        cfg.operation.times=[0,2,3,10,11,30];
        cfg.operation.factors=[ones(2,6);1,1,0,0,1,1];
        cfg.captureTimes=[0,2,2.5,3,5,8,10,11,15];
    case 'two_tuning'
        cfg.operation.times=[0,1,6,8,11,12];
        cfg.operation.factors=[ones(2,6);1,1,0.6,0.6,0.9,0.9];
        cfg.tf=12; cfg.captureTimes=[0,1,6,8,11];
    case 'two_hold'
        cfg.operation.times=[0,10]; cfg.operation.factors=ones(3,2);
        cfg.tf=10; cfg.captureTimes=[0];
    case 'single_steady'
        cfg.L=10; cfg.sectionD=0.1; cfg.nManning=0.008;
        cfg.initialDepthUp=0.095; cfg.initialDepthDown=0.095;
        cfg.QDesign=0.0013; cfg.Qofftake=0;
        cfg.prescribedFlowConstant=nan(2,1);
        cfg.operation.times=[0,10]; cfg.operation.factors=ones(2,2);
        cfg.closingDuration=[inf,inf]; cfg.finalOpening=[0,0];
        cfg.tf=10; cfg.captureTimes=0;
    case 'engineering_original'

        cfg.geometry='citygate'; cfg.L=[9000,8000,23000];
        cfg.sectionD=[8.3,8.3,8.4]; cfg.sectionH=[5.96,5.96,5.88];
        cfg.nManning=[0.015,0.015,0.015]; cfg.bedSlope=[1/4300,1/4300,1/4200];
        cfg.dx=100; cfg.initialDepthUp=[7.21,7.21,7.028];
        cfg.initialDepthDown=[7.71,7.71,7.528]; cfg.QDesign=[135,135,135];
        cfg.initialMode='original_linear_not_steady';
        cfg.HupExternal=9.5; cfg.HdownExternal=6;
        cfg.Qofftake=[0,0,0]; cfg.closingDuration=[7200,36000,36000,600];
        cfg.finalOpening=[0,0,0,0]; cfg.operation=[]; cfg.tf=43200;
        cfg.eventTimes=[600,7200,36000]; cfg.captureTimes=[0,600,3600,7200];
        cfg.prescribedFlowConstant=nan(4,1); cfg.headEps=0.02;
        cfg.storeFullField=false; cfg.progressEvery=500;
        cfg.solverOptions.flowScale=135;
        cfg.solverOptions.tolAbs=0.01; cfg.solverOptions.tolRel=1e-5;
        cfg.solverOptions.picardMaxFlowStep=30;
    otherwise
        error('IFVM:Case','Unknown caseID: %s',caseID);
end
if ~isempty(cfg.operation), cfg.eventTimes=cfg.operation.times; end
if strcmp(cfg.initialMode,'uniform_equilibrium')
    cfg.bedSlope=zeros(size(cfg.L));
    for p=1:numel(cfg.L)
        xs=CircularPipe(cfg.sectionD(p),cfg.slotRatio);
        cfg.bedSlope(p)=xs.getSf(cfg.nManning(p),cfg.initialDepthUp(p),cfg.QDesign(p));
    end
end
switch profile
    case 'screen'
    case 'paper'
        if strcmp(cfg.geometry,'circular'), cfg.dx=0.1; end
    case 'smoke'
        if strcmp(cfg.geometry,'circular'), cfg.dx=2; end
        cfg.tf=0.15; cfg.progressEvery=0; cfg.captureTimes=0;
    otherwise
        error('IFVM:Profile','Use screen, paper or smoke.');
end
if strcmp(cfg.initialMode,'uniform_equilibrium')
    cfg.HupExternal=cfg.initialDepthUp(1)-0.5*cfg.dx*cfg.bedSlope(1)+0.04;
    cfg.HdownExternal=cfg.initialDepthDown(end)+0.5*cfg.dx*cfg.bedSlope(end)-0.04;
end
end
