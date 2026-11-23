function u=CONFIG_OPENING()

u.L=[10,10];u.D=[.1,.1];u.slope=[.027,.027];u.n=[.008,.008];
u.Q0=.0013;u.dx=.1;u.slotRatio=.01;u.Cr=.5;

u.inletDepth=.13;u.initialGateHeadLoss=[.04,.04,.04];u.Cd=.83;

u.endGateStart=5;u.endGateDuration=[3,3];
u.middleGateStart=5;u.middleGateDuration=45;
u.middleGateFinalFactor=0;
u.middleTimes=[0 5 25 35 50];
u.middleFactors=[1 1 .5 .75 .5];
u.tf=60;

u.steadyChunk=5;u.steadyMaxTime=1000;
u.steadyHeadRateTol=1e-8;u.steadyFlowRateTol=1e-10;

u.holdHeadTol=1e-5;u.holdFlowTol=2e-6;
u.signFlowThreshold=1e-6;u.minimumAlternatingSteps=6;
u.tailWindow=5;u.tailHeadRangeTol=1e-5;
u.figureVisible='off';u.makeFigures=false;
u.progressEvery=2000;
u.steadyCr=.5;
u.solver.tolAbs=1e-6;u.solver.tolRel=1e-4;
u.solver.picardMaxIter=500;
u.solver.picardFixedOmega=.02;
u.solver.picardOmegaInitial=.02;
u.solver.picardOmegaMin=1e-4;u.solver.picardOmegaMax=.5;
u.headEps=2e-4;
end
