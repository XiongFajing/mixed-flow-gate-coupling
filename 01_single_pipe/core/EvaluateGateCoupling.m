function state=EvaluateGateCoupling(qGate,gateOpening,Aold,Qold,yold, ...
    nManning,bedSlope,dx,dt,xsect,Qofftake,Cd,HuExternal,HdExternal, ...
    headEps,allowReverse,context)

qGate=qGate(:); nPool=numel(yold); nGate=nPool+1;
state.valid=false; state.message=''; state.qGate=qGate;
state.A=Aold; state.Q=Qold; state.y=yold;
state.qFormula=nan(nGate,1); state.qFormulaOld=context.qOld;
state.qTarget=nan(nGate,1); state.residual=nan(nGate,1);
state.Hu=nan(nGate,1); state.Hd=nan(nGate,1);
try
    for p=1:nPool
        [M,rhs]=ApplyTrialBoundaryFlows(context.trialCache{p},qGate(p), ...
            qGate(p+1)+Qofftake(p),dx,dt);
        if isfield(context,'boundaryMomentumDepth')
            bdepth=[nan,nan];
            if p==1, bdepth(1)=context.boundaryMomentumDepth(1); end
            if p==nPool, bdepth(2)=context.boundaryMomentumDepth(2); end
            [M,rhs]=PatchBoundaryMomentum(M,rhs,yold{p},qGate(p), ...
                qGate(p+1)+Qofftake(p),xsect(p),dx,dt,bdepth,Qold{p},context.numerics.boundary);
        end
        dU=M\rhs;
        anew=Aold{p}+dU(1:2:end).'; qnew=Qold{p}+dU(2:2:end).';
        if any(~isfinite(anew)) || any(anew<=0) || any(~isfinite(qnew))
            state.message=sprintf('Pool %d produced a nonphysical state.',p); return;
        end
        ynew=xsect(p).getYofA(anew);
        if any(~isfinite(ynew))
            state.message=sprintf('Pool %d produced nonfinite heads.',p); return;
        end
        state.A{p}=anew; state.Q{p}=qnew; state.y{p}=ynew;
    end
    [state.qFormula,state.Hu,state.Hd]=GateFormulaFromState( ...
        state.y,gateOpening,HuExternal,HdExternal,bedSlope,dx,xsect, ...
        Cd,headEps,allowReverse,context.gateLaw,state.Q,context.gatePhysics);
    state.qTarget=state.qFormula;
    if isfield(context,'fixedTarget')
        fixed=isfinite(context.fixedTarget);
        state.qTarget(fixed)=context.fixedTarget(fixed);
        state.qFormula(fixed)=context.fixedInstant(fixed);
    end
    state.residual=qGate-state.qTarget;
    state.valid=all(isfinite(state.residual));
catch exception
    state.message=exception.message;
end
end
