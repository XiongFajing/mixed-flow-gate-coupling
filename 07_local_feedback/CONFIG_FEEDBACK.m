function [u,b]=CONFIG_FEEDBACK()

u=ImplicitBaseConfig();
u.Cr=.5;u.middleTimes=[0 5 50];u.middleFactors=[1 1 0];u.tf=60;
u.makeFigures=false;u.progressEvery=1000;

b.targetTime=20;
b.timeStepFactors=[1 .5 .25];
b.perturbationRatios=.01;
b.sensitivityThreshold=.1;
b.rebuildTrajectory=false;
b.csvEncoding='GB18030';
b.outputRoot=fullfile(fileparts(mfilename('fullpath')),'results');
end
