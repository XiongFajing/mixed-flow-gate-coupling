function RUN_OPENING()
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'engine'),fullfile(root,'engine','core'),'-begin');
u=CONFIG_OPENING();prep=u;prep.Cr=u.steadyCr;
c=BuildCase05(prep);[c,~]=PrepareCase05Steady(c,prep);
methods={'E1','S2'};labels={'E','I3'};
for k=1:2
 cc=c;cc.method=methods{k};cc.Cr=u.Cr;cc.dt=[];cc.headEps=u.headEps;
 fields=fieldnames(u.solver);
 for j=1:numel(fields),cc.solverOptions.(fields{j})=u.solver.(fields{j});end
 cc.quiet=true;
 timer=tic;r=RunIFVM(cc);formalSeconds=toc(timer);
 SaveResults(r,fullfile(root,'results',labels{k}),formalSeconds);
 fprintf('%s: %s, t=%.9g s, %.3f s\n',labels{k},r.status,r.t(end),formalSeconds);
end
end
