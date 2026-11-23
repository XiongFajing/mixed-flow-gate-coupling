function Qgate = OriginalGateFlow(Hup,Hdown,opening,Cd,xsect,headEps,allowReverse)

if opening<=0
    Qgate=0;
    return;
end
gateArea=min(xsect.getAofY(opening),xsect.Af);
dH=Hup-Hdown;
if allowReverse
    headFunction=sign(dH)*sqrt(abs(dH));
else
    headFunction=sqrt(max(dH,0));
end
Qgate=Cd*gateArea*sqrt(2*9.8)*headFunction;
end
