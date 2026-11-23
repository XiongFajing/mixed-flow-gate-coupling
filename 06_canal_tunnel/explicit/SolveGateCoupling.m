function [state, stats] = SolveGateCoupling( ...
    qInitial, gateOpening, Aold, Qold, yold, ...
    nManning, bedSlope, dx, dt, xsect, Qofftake, Cd, ...
    HupExternal, HdownExternal, headEps, allowReverse, options)

gateOpening = gateOpening(:);
q = qInitial(:);
closed = gateOpening <= options.closedOpeningTol;
active = find(~closed);

q(closed) = 0;
if ~allowReverse
    q(active) = max(q(active), 0);
end

evaluationCount = 0;
bestNorm = inf;
bestState = [];
bestQ = q;

previousCorrection = [];
previousNorm = inf;
omegaPreviousUpdate = options.picardOmegaInitial;
picardIteration = 0;

for picardIter = 1:options.picardMaxIter
    picardIteration = picardIter;
    current = evaluate(q);
    if ~current.valid
        break;
    end

    residualNorm = activeResidualNorm(current);
    if residualNorm < bestNorm
        bestNorm = residualNorm;
        bestState = current;
        bestQ = q;
    end

    if isConverged(current)
        state = current;
        stats = makeStats(1, 'Aitken-Picard', ...
            picardIter, 0, residualNorm, evaluationCount);
        return;
    end

    correction = -current.residual(active);

    if isempty(previousCorrection)
        omega = options.picardOmegaInitial;
    else
        deltaCorrection = correction - previousCorrection;
        denominator = deltaCorrection.'*deltaCorrection;

        if denominator > 1e-24
            omegaCandidate = -omegaPreviousUpdate ...
                *(previousCorrection.'*deltaCorrection)/denominator;
        else
            omegaCandidate = 0.5*omegaPreviousUpdate;
        end

        if ~isfinite(omegaCandidate) || omegaCandidate <= 0
            omegaCandidate = 0.5*omegaPreviousUpdate;
        end

        omega = min(max(omegaCandidate, ...
            options.picardOmegaMin), options.picardOmegaMax);

        if residualNorm > 1.2*previousNorm
            omega = max(options.picardOmegaMin, ...
                min(omega, 0.5*omegaPreviousUpdate));
        end
    end

    flowIncrement = omega*correction;
    stepNorm = norm(flowIncrement, inf);
    if stepNorm > options.picardMaxFlowStep
        flowIncrement = flowIncrement ...
            *options.picardMaxFlowStep/stepNorm;
    end

    q(active) = q(active) + flowIncrement;
    q(closed) = 0;
    if ~allowReverse
        q(active) = max(q(active), 0);
    end

    previousCorrection = correction;
    previousNorm = residualNorm;
    omegaPreviousUpdate = omega;
end

if isempty(bestState)
    initialState = evaluate(q);
    if ~initialState.valid
        error(['Gate coupling could not evaluate a physical trial state. ', ...
            'The previous time state was not overwritten. Detail: %s'], ...
            initialState.message);
    end
    bestState = initialState;
    bestQ = q;
    bestNorm = activeResidualNorm(initialState);
end

q = bestQ;
newtonIteration = 0;

for newtonIter = 1:options.newtonMaxIter
    newtonIteration = newtonIter;
    base = evaluate(q);
    if ~base.valid
        break;
    end

    baseNorm = activeResidualNorm(base);
    if baseNorm < bestNorm
        bestNorm = baseNorm;
        bestState = base;
        bestQ = q;
    end

    if isConverged(base)
        state = base;
        stats = makeStats(2, 'Damped-Newton', ...
            picardIteration, newtonIter, baseNorm, evaluationCount);
        return;
    end

    nActive = numel(active);
    Jacobian = zeros(nActive, nActive);
    baseResidual = base.residual(active);
    jacobianValid = true;

    for column = 1:nActive
        gate = active(column);
        perturbation = options.newtonFdAbs ...
            + options.newtonFdRel*max(abs(q(gate)), 1);

        qPerturbed = q;
        qPerturbed(gate) = qPerturbed(gate) + perturbation;
        perturbed = evaluate(qPerturbed);

        if perturbed.valid
            Jacobian(:, column) = ...
                (perturbed.residual(active) - baseResidual)/perturbation;
        else

            qPerturbed = q;
            qPerturbed(gate) = qPerturbed(gate) - perturbation;
            if ~allowReverse
                qPerturbed(gate) = max(qPerturbed(gate), 0);
            end
            actualPerturbation = qPerturbed(gate) - q(gate);
            perturbed = evaluate(qPerturbed);

            if ~perturbed.valid || abs(actualPerturbation) < 1e-14
                jacobianValid = false;
                break;
            end
            Jacobian(:, column) = ...
                (perturbed.residual(active) - baseResidual) ...
                /actualPerturbation;
        end
    end

    if ~jacobianValid || any(~isfinite(Jacobian(:)))
        break;
    end

    if rcond(Jacobian) > 1e-11
        newtonStep = -Jacobian\baseResidual;
    else
        regularization = 1e-8*max(1, norm(Jacobian, 'fro')^2);
        newtonStep = -(Jacobian.'*Jacobian ...
            + regularization*eye(nActive))\(Jacobian.'*baseResidual);
    end

    stepNorm = norm(newtonStep, inf);
    if stepNorm > options.newtonMaxFlowStep
        newtonStep = newtonStep ...
            *options.newtonMaxFlowStep/stepNorm;
    end

    accepted = false;
    for lineSearch = 0:options.lineSearchMaxIter
        alpha = 0.5^lineSearch;
        qCandidate = q;
        qCandidate(active) = q(active) + alpha*newtonStep;
        qCandidate(closed) = 0;
        if ~allowReverse
            qCandidate(active) = max(qCandidate(active), 0);
        end

        candidate = evaluate(qCandidate);
        if ~candidate.valid
            continue;
        end

        candidateNorm = activeResidualNorm(candidate);
        if candidateNorm < baseNorm*(1 - 1e-4*alpha) ...
                || isConverged(candidate)
            q = qCandidate;
            accepted = true;

            if candidateNorm < bestNorm
                bestNorm = candidateNorm;
                bestState = candidate;
                bestQ = qCandidate;
            end
            break;
        end
    end

    if ~accepted
        break;
    end
end

if ~isempty(bestState) && isConverged(bestState)
    state = bestState;
    stats = makeStats(2, 'Damped-Newton', ...
        picardIteration, newtonIteration, bestNorm, evaluationCount);
    return;
end

[worstResidual, localIndex] = max(abs(bestState.residual(active)));
worstGate = active(localIndex);
error(['Gate coupling failed: max residual %.6g m^3/s at gate %d. ', ...
    'The previous time state was NOT overwritten. Suggested checks: ', ...
    'halve dt, increase newtonMaxIter, or inspect the gate/head data.'], ...
    worstResidual, worstGate);

    function evaluated = evaluate(qTrial)
        evaluated = EvaluateGateCoupling( ...
            qTrial, gateOpening, Aold, Qold, yold, ...
            nManning, bedSlope, dx, dt, xsect, Qofftake, Cd, ...
            HupExternal, HdownExternal, headEps, allowReverse);
        evaluationCount = evaluationCount + 1;
    end

    function value = activeResidualNorm(evaluated)
        if isempty(active)
            value = 0;
        else
            value = norm(evaluated.residual(active), inf);
        end
    end

    function yes = isConverged(evaluated)
        if ~evaluated.valid
            yes = false;
            return;
        end
        if isempty(active)
            yes = true;
            return;
        end
        scale = max(1, norm(evaluated.qFormula(active), inf));
        tolerance = options.tolAbs + options.tolRel*scale;
        yes = activeResidualNorm(evaluated) <= tolerance;
    end

    function output = makeStats( ...
            methodCode, methodName, picardCount, newtonCount, ...
            residualNorm, evaluations)
        output.methodCode = methodCode;
        output.methodName = methodName;
        output.picardIterations = picardCount;
        output.newtonIterations = newtonCount;
        output.iterations = picardCount + newtonCount;
        output.residualNorm = residualNorm;
        output.evaluations = evaluations;
    end
end
