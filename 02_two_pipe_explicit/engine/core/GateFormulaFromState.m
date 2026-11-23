function [qInstant, Hu, Hd] = GateFormulaFromState( ...
    y, gateOpening, HupExternal, HdownExternal, ...
    bedSlope, dx, xsect, Cd, headEps, allowReverse, gateLaw, Qcell, physics)

if nargin < 11, gateLaw = 'smooth'; end
if nargin < 12,Qcell=[];end
if nargin < 13,physics='submerged';end
nPool = numel(y);
nGate = nPool + 1;
gateOpening = gateOpening(:);
if numel(gateOpening) ~= nGate
    error('IFVM:GateSize', 'Gate opening must have nPool+1 entries.');
end
[Hu, Hd] = GetGateHeads(y, HupExternal, HdownExternal, bedSlope, dx);
qInstant = zeros(nGate, 1);
for gate = 1:nGate
    sectionIndex = max(1, min(gate - 1, nPool));
    if strcmp(physics,'free_outlet') && gate==nGate && gateOpening(gate)>0
        if gateOpening(gate)>=xsect(sectionIndex).D && ~isempty(Qcell)

            qInstant(gate)=Qcell{end}(end);
            if Hd(gate)<=0,qInstant(gate)=max(qInstant(gate),0);end
            continue;
        end
        Hd(gate)=max(Hd(gate),0);
        effectiveTail=max(Hd(gate),0.5*gateOpening(gate));
        if Hd(gate)<=0 && Hu(gate)<=gateOpening(gate) && ~isempty(Qcell)
            qInstant(gate)=max(Qcell{end}(end),0);continue;
        end
        reverseHere=allowReverse && Hd(gate)>0;
        qInstant(gate)=GateFlow(Hu(gate),effectiveTail,gateOpening(gate), ...
            Cd,xsect(sectionIndex),headEps,reverseHere,gateLaw);
        continue;
    end
    qInstant(gate) = GateFlow( ...
        Hu(gate), Hd(gate), gateOpening(gate), ...
        Cd, xsect(sectionIndex), headEps, allowReverse, gateLaw);
end
end
