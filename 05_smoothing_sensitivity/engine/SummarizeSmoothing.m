function SummarizeSmoothing(folder)
s=load(fullfile(folder,'batch_config.mat'),'u','b');u=s.u;b=s.b;
s=load(fullfile(folder,'initial_state.mat'),'Href','Qref');Href=s.Href;Qref=s.Qref;
refs=cell(1,2);refOK=false(1,2);
for j=1:2
 f=fullfile(folder,sprintf('reference_%d',j),'results.mat');
 if isfile(f)
  s=load(f,'r','meta');refs{j}=s.r;
  refOK(j)=strcmp(s.meta.status,'completed') && ~isempty(s.r);
 end
end
refQ=NaN;refH=NaN;valid=false;reason='参考未完成或未启用';
if b.referenceEnabled && all(refOK)
 [eq,eh]=Differences(refs{1},refs{2},b.errorWindow);
 refQ=100*eq/Qref;refH=100*eh/Href;
 valid=refQ<=b.referenceAgreementPercent && refH<=b.referenceAgreementPercent;
 if valid,reason='两级参考一致性通过';else,reason='两级参考差异超限，先不要使用误差曲线';end
end
T=table(refOK(1),refOK(2),refQ,refH,b.referenceAgreementPercent,valid,string(reason), ...
 'VariableNames',{'Reference1Completed','Reference2Completed','FlowDifference_percent','HeadDifference_percent', ...
 'Tolerance_percent','ReferenceValid','Message'});
writetable(T,fullfile(folder,'reference_check.csv'),'Encoding','UTF-8');
n=numel(b.multipliers);
M=nan(n,25);status=strings(n,1);message=strings(n,1);
for j=1:n
 epsH=b.baseHeadEps*b.multipliers(j);
 M(j,1:4)=[j,b.multipliers(j),epsH,epsH/Href];status(j)="未运行";
 f=fullfile(folder,sprintf('group_%02d',j),'results.mat');
 if ~isfile(f),continue;end
 s=load(f,'r','meta');r=s.r;meta=s.meta;
 status(j)=string(meta.status);message(j)=string(meta.message);
 M(j,10:11)=[meta.formalSeconds,meta.exportSeconds];
 if isempty(r),continue;end
 M(j,5:9)=[r.t(end),u.tf,100*r.t(end)/u.tf,r.completedSteps,r.plannedSteps];
 if ~isempty(r.iterations)
  M(j,12:15)=[mean(r.iterations),max(r.iterations),sum(r.iterations),r.totalEvaluations];
 else
  M(j,14:15)=[0,r.totalEvaluations];
 end
 M(j,16)=r.computeSeconds;
 M(j,17)=r.initialHeadCompatibilityError;
 M(j,18)=max([0;abs(r.massResidual(:))]);
 if isfield(r.failure,'timeStart'),M(j,19:20)=[r.failure.timeStart,r.failure.timeEnd];end
 if strcmp(r.status,'completed') && valid
  [eq,eh,np]=Differences(r,refs{2},b.errorWindow);
  M(j,21:25)=[eq,eh,100*eq/Qref,100*eh/Href,np];
 end
end
names={'Group','Multiplier','HeadEps_m','Alpha', ...
 'EndTime_s','PlannedTime_s','Completion_percent','CompletedSteps','PlannedSteps', ...
 'FormalTime_s','ExportTime_s','MeanEvaluations','MaxEvaluations', ...
 'TotalAcceptedEvaluations','TotalEvaluations','LoopTime_s', ...
 'InitialFlowMismatch_m3s','MaxMassResidual_m3','FailureStart_s','FailureEnd_s', ...
 'Gate2_FlowRMSE_m3s','Gate2_HeadDifferenceRMSE_m','Gate2_FlowNRMSE_percent','Gate2_HeadDifferenceNRMSE_percent','ErrorPoints'};
T = array2table(M,'VariableNames',names);
T.('Status') = status;
T.('Message') = message;
T=T(:,{'Group','HeadEps_m','Alpha','Status','EndTime_s','FailureStart_s','MeanEvaluations','MaxEvaluations','FormalTime_s','Gate2_FlowNRMSE_percent','Gate2_HeadDifferenceNRMSE_percent'});
writetable(T,fullfile(folder,'metrics.csv'),'Encoding','UTF-8');
save(fullfile(folder,'metrics.mat'),'T','Href','Qref','valid','refQ','refH');
end

function [eq,eh,np]=Differences(r,ref,window)

t=r.t(2:end);tr=ref.t(2:end);
assert(numel(t)==numel(tr) && max(abs(t-tr))<1e-10, ...
 '结果与参考时间网格不一致，不能直接比较。');
mask=t>=window(1)-1e-10 & t<=window(2)+1e-10;
assert(any(mask),'误差统计时段没有接受步末数据。');
dq=r.qBoundary(2,mask)-ref.qBoundary(2,mask);
h=r.gateHu(2,2:end)-r.gateHd(2,2:end);
hr=ref.gateHu(2,2:end)-ref.gateHd(2,2:end);
dh=h(mask)-hr(mask);np=nnz(mask);
eq=sqrt(mean(dq.^2));eh=sqrt(mean(dh.^2));
end
