function RUN_IMPLICIT()
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'engine'),fullfile(root,'engine','core'),'-begin');
u=CONFIG_IMPLICIT();prep=u;prep.Cr=u.steadyCr;
c=BuildCase05(prep);[c,~]=PrepareCase05Steady(c,prep);
methods={'R1','S1','S2'};labels={'I1','I2','I3'};
for k=1:3
 cc=c;cc.method=methods{k};cc.Cr=u.Cr;cc.dt=[];cc.headEps=u.headEps;
 f=fieldnames(u.solver);
 for j=1:numel(f),cc.solverOptions.(f{j})=u.solver.(f{j});end
 cc.quiet=true;
 timer=tic;r=RunIFVM(cc);formalSeconds=toc(timer);
 SaveResults(r,fullfile(root,'results',labels{k}),formalSeconds);
 fprintf('%s: %s, t=%.9g s, %.3f s\n',labels{k},r.status,r.t(end),formalSeconds);
end
end
