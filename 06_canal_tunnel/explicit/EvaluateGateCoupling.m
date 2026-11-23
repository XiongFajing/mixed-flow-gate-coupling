function state = EvaluateGateCoupling( ...
    qGate, gateOpening, Aold, Qold, yold, ...
    nManning, bedSlope, dx, dt, xsect, Qofftake, Cd, ...
    HupExternal, HdownExternal, headEps, allowReverse)

nPool = numel(yold);
nGate = nPool + 1;
qGate = qGate(:);
gateOpening = gateOpening(:);

state.valid = false;
state.message = '';
state.qGate = qGate;
state.qFormula = nan(nGate, 1);
state.residual = nan(nGate, 1);
state.A = Aold;
state.Q = Qold;
state.y = yold;
state.Hu = nan(nGate, 1);
state.Hd = nan(nGate, 1);

if numel(qGate) ~= nGate || numel(gateOpening) ~= nGate
    state.message = 'Gate vector has an inconsistent size.';
    return;
end

try
    for pool = 1:nPool
        QinBoundary = qGate(pool);
        QoutBoundary = qGate(pool + 1) + Qofftake(pool);

        [Mmat, rhs] = FVMMatrixNoGhost( ...
            yold{pool}, Qold{pool}, nManning(pool), bedSlope(pool), ...
            dx, dt, xsect(pool), QinBoundary, QoutBoundary);

        dU = Mmat\rhs;
        Anew = Aold{pool} + dU(1:2:end).';
        Qnew = Qold{pool} + dU(2:2:end).';

        if any(~isfinite(Anew)) || any(Anew <= 0) || any(~isfinite(Qnew))
            state.message = sprintf('Pool %d produced a nonphysical state.', pool);
            return;
        end

        ynew = xsect(pool).getYofA(Anew);
        if any(~isfinite(ynew))
            state.message = sprintf('Pool %d produced a non-finite head.', pool);
            return;
        end

        state.A{pool} = Anew;
        state.Q{pool} = Qnew;
        state.y{pool} = ynew;
    end

    [state.Hu, state.Hd] = GetGateHeads( ...
        state.y, HupExternal, HdownExternal, bedSlope, dx);

    for gate = 1:nGate
        if gate == 1
            sectionIndex = 1;
        elseif gate == nGate
            sectionIndex = nPool;
        else
            sectionIndex = gate - 1;
        end

        state.qFormula(gate) = OriginalGateFlow( ...
            state.Hu(gate), state.Hd(gate), gateOpening(gate), ...
            Cd, xsect(sectionIndex), headEps, allowReverse);
    end

    state.residual = qGate - state.qFormula;
    state.valid = all(isfinite(state.residual));
    if ~state.valid
        state.message = 'Gate residual contains a non-finite value.';
    end
catch exception
    state.message = exception.message;
end
end
