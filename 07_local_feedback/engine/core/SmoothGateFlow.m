function Qgate = SmoothGateFlow( ...
    Hup, Hdown, opening, Cd, xsect, headEps, allowReverse)

if opening <= 0
    Qgate = 0;
    return;
end

gateArea = min(xsect.getAofY(opening), xsect.Af);
dH = Hup - Hdown;
denominator = (dH^2 + headEps^2)^(1/4);

if allowReverse
    headFunction = dH/denominator;
else
    headFunction = max(dH, 0)/denominator;
end

Qgate = Cd*gateArea*sqrt(2*9.8)*headFunction;
end
