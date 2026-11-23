function [Flux, Jacobian] = BoundaryFlux(yCell, QBoundary, xsect,QCell,side,kind)

A = xsect.getAofY(yCell);
I1 = xsect.getIofY(yCell);
c = xsect.getCofY(yCell);
uBoundary = QBoundary/A;

Flux = [ ...
    QBoundary; ...
    QBoundary^2/A + 9.8*I1];

Jacobian = [ ...
    0, 0; ...
    c^2 - uBoundary^2, 0];
if nargin>=6 && strcmp(kind,'riemann')

    sigma=max(abs(QCell),abs(2*QBoundary-QCell))/A+c;
    V=QCell^2-2*QBoundary*QCell+2*QBoundary^2;
    Flux(2)=V/A+9.8*I1+side*sigma*(QCell-QBoundary);
    Jacobian(2,1)=c^2-V/A^2;
    Jacobian(2,2)=2*(QCell-QBoundary)/A+side*sigma;
end
end
