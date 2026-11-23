function [u,b]=CONFIG_SMOOTHING()

u=ImplicitBaseConfig();

u.middleTimes=[0 5 50]; u.middleFactors=[1 1 0];

u.tf=60; u.makeFigures=false; u.progressEvery=2000;
b.method='S2';
b.baseHeadEps=2e-4;
b.multipliers=[.01 .03 .1 .3 1 3 10 30 100];

b.referenceEnabled=true;
b.referenceHeadEps=[2e-8 1e-8];
b.referenceTolAbs=1e-10; b.referenceTolRel=1e-8;
b.referenceMaxIter=2500;
b.referenceAgreementPercent=.1;
b.referenceMinEpsRatio=10;
b.errorWindow=[0 u.tf];
b.resume=true;
b.rerunFailed=false;
b.forceSteady=false;
b.outputRoot=fullfile(fileparts(mfilename('fullpath')),'results');

end
