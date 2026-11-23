function ExportEngineering(folder,methodName,makeFigures,ts,hu,hd,qu,qd,t0,t1,G,q,qoff,seconds,evals,iters,newtons,residuals)
data=ts(:);names={'Time_s'};
for p=1:size(hu,1)
 data=[data,hu(p,:).',hd(p,:).',qu(p,:).',qd(p,:).'];
 names=[names,{sprintf('Pool%d_UpstreamHead_m',p),sprintf('Pool%d_DownstreamHead_m',p),sprintf('Pool%d_UpstreamCellFlow_m3s',p),sprintf('Pool%d_DownstreamCellFlow_m3s',p)}];
end
WriteNumeric(fullfile(folder,'states.csv'),names,data);
data=[t0(:),t1(:)];names={'StartTime_s','EndTime_s'};
for g=1:size(G,1)
 data=[data,G(g,:).',q(g,:).'];names=[names,{sprintf('Gate%d_Opening_m',g),sprintf('Gate%d_Flow_m3s',g)}];
end
WriteNumeric(fullfile(folder,'gate_flows.csv'),names,data);
save(fullfile(folder,'results.mat'),'ts','hu','hd','qu','qd','t0','t1','G','q','qoff','evals','iters','newtons','residuals','seconds');
end

function WriteNumeric(file,names,data)
fid=fopen(file,'w','n','UTF-8');
assert(fid>=0,'Cannot write file: %s',file);
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,uint8([239 187 191]),'uint8');
fprintf(fid,'%s\n',strjoin(names,','));
if ~isempty(data)
 fmt=[repmat('%.16g,',1,size(data,2)-1),'%.16g\n'];
 fprintf(fid,fmt,data.');
end
end
