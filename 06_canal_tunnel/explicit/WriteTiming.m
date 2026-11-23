function WriteTiming(folder,method,status,message,Cr,dx,dt,tf,t,n,init,loop,solve,output,evals,iters,newtonSteps)
T=table(string(method),string(status),Cr,dx,dt,tf,t,n,init,solve,loop,output,evals,iters,newtonSteps,string(message), ...
 'VariableNames',{'Method','Status','Cr','Dx_m','NominalDt_s','PlannedTime_s','EndTime_s','CompletedSteps','InitializationTime_s','SolverTime_s','LoopTime_s','ExportTime_s','TotalEvaluations','TotalIterations','NewtonSteps','Message'});
writetable(T,fullfile(folder,'summary.csv'));
end
