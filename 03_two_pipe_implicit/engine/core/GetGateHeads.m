function [Hu, Hd] = GetGateHeads( ...
    yCell, HupExternal, HdownExternal, bedSlope, dx)

nPool = numel(yCell);
Hu = zeros(nPool + 1, 1);
Hd = zeros(nPool + 1, 1);

Hu(1) = HupExternal;
Hd(1) = yCell{1}(1) - 0.5*dx*bedSlope(1);

for gate = 2:nPool

    Hu(gate) = yCell{gate - 1}(end) ...
        + 0.5*dx*bedSlope(gate - 1);

    Hd(gate) = yCell{gate}(1) ...
        - 0.5*dx*bedSlope(gate);
end

Hu(nPool + 1) = yCell{nPool}(end) ...
    + 0.5*dx*bedSlope(nPool);
Hd(nPool + 1) = HdownExternal;
end
