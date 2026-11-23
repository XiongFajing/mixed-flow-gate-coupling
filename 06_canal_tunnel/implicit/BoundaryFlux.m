function [Flux, Jacobian] = BoundaryFlux(yCell, QBoundary, xsect)

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
end
