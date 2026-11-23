function [state,stats]=SolveGateCoupling(qInitial,opening,Aold,Qold,yold, ...
    nManning,bedSlope,dx,dt,xsect,Qofftake,Cd,HuExternal,HdExternal, ...
    headEps,allowReverse,options)

settings=MethodSettings(options.method);
opening=opening(:); q=qInitial(:); nGate=numel(yold)+1;
assert(numel(opening)==nGate && numel(q)==nGate,'Gate vector size mismatch.');
assert(options.picardMaxIter>=1 && options.picardMaxIter==floor(options.picardMaxIter));
assert(options.picardFixedOmega>0 && options.picardFixedOmega<=1);
context.gateLaw=settings.gateLaw;
if ~isfield(options,'gatePhysics'),options.gatePhysics='submerged';end
context.gatePhysics=options.gatePhysics;
context.timeCoupling=settings.timeCoupling;
context.qOld=GateFormulaFromState(yold,options.gateOpeningOld, ...
    options.externalHeadOld(1),options.externalHeadOld(2),bedSlope,dx,xsect,Cd,headEps,allowReverse,settings.gateLaw,Qold,context.gatePhysics);
context.fixedTarget=nan(nGate,1);
if ~isfield(options,'fixedFlowOld'), options.fixedFlowOld=nan(nGate,1); options.fixedFlowNew=nan(nGate,1); end
fixed=isfinite(options.fixedFlowOld);
if ~isequal(fixed,isfinite(options.fixedFlowNew)), error('IFVM:BoundaryType','Prescribed boundary type changes within a step.'); end
context.qOld(fixed)=options.fixedFlowOld(fixed);
context.fixedTarget(fixed)=options.fixedFlowNew(fixed);
context.fixedInstant=options.fixedFlowNew;

context.boundaryMomentumDepth=options.boundaryMomentumDepth;
if ~isfield(options,'numerics'),options.numerics.kernel='rusanov';end
if ~isfield(options.numerics,'boundary'),options.numerics.boundary='direct';end
context.numerics=options.numerics;
context.trialCache=BuildIFVMTrialCache(yold,Qold,nManning,bedSlope,dx,dt,xsect,options.numerics);
closed=opening<=options.closedOpeningTol; active=find(~closed & ~fixed);
stats.methodName=settings.method; stats.methodCode=settings.number; stats.iterations=0;
stats.evaluations=0; stats.residualNorm=inf; stats.converged=false;
stats.message=''; stats.lastQ=q; stats.lastResidual=nan(nGate,1);
stats.residualHistory=[]; stats.omegaHistory=[];

stats.iterationTrace=zeros(0,5);
state=[];
recordConvergence=isfield(options,'recordConvergence') && options.recordConvergence;
if recordConvergence
    stats.convergenceTrace=struct();
end
if strcmp(settings.relaxation,'none')

    q=context.qOld;
    if strcmp(settings.timeCoupling,'explicit_old_new_opening')
        q=GateFormulaFromState(yold,opening,options.externalHeadOld(1), ...
            options.externalHeadOld(2),bedSlope,dx,xsect,Cd,headEps,allowReverse,settings.gateLaw,Qold,context.gatePhysics);
        q(fixed)=options.fixedFlowNew(fixed);
    end
    trial=evaluate(q);
    stats.lastQ=q; stats.lastResidual=trial.residual;
    if trial.valid
        state=trial; stats.converged=true;
        stats.residualNorm=norm(trial.residual,inf);
    else
        stats.message=trial.message;
    end
    return;
end
q(closed)=0; q(fixed)=context.fixedTarget(fixed);
if ~allowReverse, q(active)=max(q(active),0); end
previousCorrection=[]; previousNorm=inf;
omegaPrevious=options.picardOmegaInitial;
for m=1:options.picardMaxIter
    stats.iterations=m;
    trial=evaluate(q);
    stats.lastQ=q; stats.lastResidual=trial.residual;
    if ~trial.valid
        stats.message=trial.message; return;
    end
    if isempty(active), residualNorm=0; else, residualNorm=norm(trial.residual(active),inf); end
    stats.residualNorm=residualNorm;
    stats.residualHistory(m)=residualNorm;
    stats.lastHu=trial.Hu; stats.lastHd=trial.Hd;
    if isempty(active), scale=options.flowScale;
    else, scale=max(options.flowScale,norm(trial.qTarget(active),inf)); end
    if recordConvergence
        stats.convergenceTrace.iteration(m,1)=m;
        stats.convergenceTrace.residual(m,1)=residualNorm;
        stats.convergenceTrace.tolerance(m,1)=options.tolAbs+options.tolRel*scale;
        stats.convergenceTrace.normalized(m,1)=residualNorm/stats.convergenceTrace.tolerance(m,1);
        stats.convergenceTrace.Hu(m,:)=trial.Hu(:).';
        stats.convergenceTrace.Hd(m,:)=trial.Hd(:).';
        stats.convergenceTrace.dH(m,:)=(trial.Hu(:)-trial.Hd(:)).';
        stats.convergenceTrace.qTrial(m,:)=q(:).';
        stats.convergenceTrace.qTarget(m,:)=trial.qTarget(:).';
        stats.convergenceTrace.gateResidual(m,:)=trial.residual(:).';
        stats.convergenceTrace.omega(m,1)=NaN;
    end
    if residualNorm<=options.tolAbs+options.tolRel*scale
        state=trial; stats.converged=true; return;
    end
    if m==options.picardMaxIter, break; end
    correction=-trial.residual(active);
    if strcmp(settings.relaxation,'fixed')
        omega=options.picardFixedOmega;
    elseif isempty(previousCorrection)
        omega=options.picardOmegaInitial;
    else
        delta=correction-previousCorrection; den=delta.'*delta;
        if den>1e-24
            omega=-omegaPrevious*(previousCorrection.'*delta)/den;
        else
            omega=0.5*omegaPrevious;
        end
        if ~isfinite(omega) || omega<=0, omega=0.5*omegaPrevious; end
        omega=min(max(omega,options.picardOmegaMin),options.picardOmegaMax);
        if residualNorm>1.2*previousNorm
            omega=max(options.picardOmegaMin,min(omega,0.5*omegaPrevious));
        end
    end
    increment=omega*correction; stepNorm=norm(increment,inf);
    stats.omegaHistory(m)=omega;
    if recordConvergence,stats.convergenceTrace.omega(m,1)=omega;end
    if stepNorm>options.picardMaxFlowStep
        increment=increment*options.picardMaxFlowStep/stepNorm;
    end
    q(active)=q(active)+increment; q(closed)=0; q(fixed)=context.fixedTarget(fixed);
    if ~allowReverse, q(active)=max(q(active),0); end
    previousCorrection=correction; previousNorm=residualNorm; omegaPrevious=omega;
end
stats.message=sprintf('Picard did not converge in %d evaluations; residual %.9g m^3/s.', ...
    options.picardMaxIter,stats.residualNorm);
    function trial=evaluate(qTrial)
        trial=EvaluateGateCoupling(qTrial,opening,Aold,Qold,yold,nManning, ...
            bedSlope,dx,dt,xsect,Qofftake,Cd,HuExternal,HdExternal, ...
            headEps,allowReverse,context);
        stats.evaluations=stats.evaluations+1;
        if trial.valid
            if isempty(active)
                rn=0;
                scale=options.flowScale;
            else
                rn=norm(trial.residual(active),inf);
                scale=max(options.flowScale,norm(trial.qTarget(active),inf));
            end

            threshold=options.tolAbs+options.tolRel*scale;

            stats.iterationTrace(end+1,:)=[ ...
                stats.evaluations, ...
                rn, ...
                threshold, ...
                trial.Hu(2)-trial.Hd(2), ...
                qTrial(2)];
        end
    end
end
