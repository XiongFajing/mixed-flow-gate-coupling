function [state, stats] = AdvanceGateNoPicard( ...
    gateOpening, Aold, Qold, yold, ...
    nManning, bedSlope, dx, dt, xsect, Qofftake, Cd, ...
    HupExternal, HdownExternal, headEps, allowReverse)

nPool = numel(yold);
nGate = nPool + 1;
gateOpening = gateOpening(:);

if numel(gateOpening) ~= nGate
    error('Gate opening vector has an inconsistent size.');
end

[HuOld, HdOld] = GetGateHeads( ...
    yold, HupExternal, HdownExternal, bedSlope, dx);

qGate = zeros(nGate, 1);
for gate = 1:nGate
    if gate == 1
        sectionIndex = 1;
    elseif gate == nGate
        sectionIndex = nPool;
    else
        sectionIndex = gate - 1;
    end

    qGate(gate) = SmoothGateFlow( ...
        HuOld(gate), HdOld(gate), gateOpening(gate), ...
        Cd, xsect(sectionIndex), headEps, allowReverse);
end

state = EvaluateGateCoupling( ...
    qGate, gateOpening, Aold, Qold, yold, ...
    nManning, bedSlope, dx, dt, xsect, Qofftake, Cd, ...
    HupExternal, HdownExternal, headEps, allowReverse);

if ~state.valid
    error(['No-Picard step produced an invalid hydraulic state. ', ...
        'The previous time state was not overwritten. Detail: %s'], ...
        state.message);
end

state.HuOld = HuOld;
state.HdOld = HdOld;

stats.methodCode = 0;
stats.methodName = 'No-Picard explicit gate';
stats.picardIterations = 0;
stats.newtonIterations = 0;
stats.iterations = 0;
stats.evaluations = 1;
stats.residualNorm = norm(state.residual, inf);
end
