function RUN_EXPLICIT()
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'engine'),fullfile(root,'engine','core'),'-begin');
values=[.5 .25 .125];
for k=1:numel(values)
 u=CONFIG_EXPLICIT();u.Cr=values(k);
 c=BuildCase05(u);[c,~]=PrepareCase05Steady(c,u);
 c.quiet=true;
 timer=tic;r=RunIFVM(c);formalSeconds=toc(timer);
 folder=fullfile(root,'results',sprintf('Cr_%g',values(k)));
 SaveResults(r,folder,formalSeconds);
 fprintf('Cr=%g: %s, t=%.9g s, %.3f s\n',values(k),r.status,r.t(end),formalSeconds);
end
end
