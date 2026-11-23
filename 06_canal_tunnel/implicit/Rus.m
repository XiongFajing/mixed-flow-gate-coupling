function [Flux, spectralRadius] = Rus(yl, yr, Ql, Qr, xsect)

Al = xsect.getAofY(yl);
Ar = xsect.getAofY(yr);
cl = xsect.getCofY(yl);
cr = xsect.getCofY(yr);
Il = xsect.getIofY(yl);
Ir = xsect.getIofY(yr);

ul = Ql/Al;
ur = Qr/Ar;

spectralRadius = max(abs(ul) + cl, abs(ur) + cr);

Ul = [Al; Ql];
Ur = [Ar; Qr];
FluxLeft = [Ql; Ql^2/Al + 9.8*Il];
FluxRight = [Qr; Qr^2/Ar + 9.8*Ir];

Flux = 0.5*( ...
    FluxLeft + FluxRight + spectralRadius*(Ul - Ur));
end
