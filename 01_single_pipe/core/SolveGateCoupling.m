function [state,stats]=SolveGateCoupling(qInitial,opening,Aold,Qold,yold, ...
    nManning,bedSlope,dx,dt,xsect,Qofftake,Cd,HuExternal,HdExternal, ...
    headEps,allowReverse,options)

settings=MethodSettings(options.method);
opening=opening(:); q=qInitial(:); nGate=numel(yold)+1;
assert(numel(opening)==nGate && numel(q)==nGate,'Gate vector size mismatch.');
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
stats.methodName=settings.method; stats.methodCode=settings.number; stats.iterations=0;
stats.evaluations=0; stats.residualNorm=inf; stats.converged=false;
stats.message=''; stats.lastQ=q; stats.lastResidual=nan(nGate,1);
stats.residualHistory=[]; stats.omegaHistory=[];
state=[];
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
    function trial=evaluate(qTrial)
        trial=EvaluateGateCoupling(qTrial,opening,Aold,Qold,yold,nManning, ...
            bedSlope,dx,dt,xsect,Qofftake,Cd,HuExternal,HdExternal, ...
            headEps,allowReverse,context);
        stats.evaluations=stats.evaluations+1;
    end
end
