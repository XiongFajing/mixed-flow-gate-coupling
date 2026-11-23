function [Aeq, Beq] = FVMMatrixNoGhost( ...
    y, Q, nManning, S0, dx, dt, xsect, QinBoundary, QoutBoundary)

N = numel(y);
if N < 2
    error('At least two real control volumes are required.');
end

if numel(Q) ~= N
    error('y and Q must have the same length.');
end

lambda = dt/dx;
Aeq = spalloc(2*N, 2*N, 12*N);
Beq = zeros(2*N, 1);
I2 = eye(2);

for j = 1:N
    yj = y(j);
    Qj = Q(j);

    Aj = xsect.getAofY(yj);
    cj = xsect.getCofY(yj);
    uj = Qj/Aj;
    Jj = [0, 1; cj^2 - uj^2, 2*uj];

    Sf = xsect.getSf(nManning, yj, Qj);
    dSfdA = xsect.getdSfdA(nManning, yj, Qj);
    dSfdQ = xsect.getdSfdQ(nManning, yj, Qj);

    sourceJac = zeros(2, 2);
    sourceJac(2, 1) = 9.8*(S0 - Sf - Aj*dSfdA);
    sourceJac(2, 2) = -9.8*Aj*dSfdQ;

    source = [0; 9.8*Aj*(S0 - Sf)];
    row = 2*j-1:2*j;

    if j == 1
        [Fleft, Jleft] = BoundaryFlux(yj, QinBoundary, xsect);

        yr = y(j + 1);
        Qr = Q(j + 1);
        Ar = xsect.getAofY(yr);
        cr = xsect.getCofY(yr);
        ur = Qr/Ar;
        Jr = [0, 1; cr^2 - ur^2, 2*ur];

        [Fright, spectralRight] = Rus(yj, yr, Qj, Qr, xsect);

        BB = I2 ...
            - lambda*Jleft ...
            + 0.5*lambda*(Jj + spectralRight*I2) ...
            - dt*sourceJac;
        CC = 0.5*lambda*(Jr - spectralRight*I2);

        Aeq(row, 1:4) = [BB, CC];
        Beq(row) = lambda*(Fleft - Fright) + dt*source;

    elseif j == N
        yl = y(j - 1);
        Ql = Q(j - 1);
        Al = xsect.getAofY(yl);
        cl = xsect.getCofY(yl);
        ul = Ql/Al;
        Jl = [0, 1; cl^2 - ul^2, 2*ul];

        [Fleft, spectralLeft] = Rus(yl, yj, Ql, Qj, xsect);
        [Fright, Jright] = BoundaryFlux(yj, QoutBoundary, xsect);

        AA = -0.5*lambda*(Jl + spectralLeft*I2);
        BB = I2 ...
            - 0.5*lambda*(Jj - spectralLeft*I2) ...
            + lambda*Jright ...
            - dt*sourceJac;

        Aeq(row, 2*N-3:2*N) = [AA, BB];
        Beq(row) = lambda*(Fleft - Fright) + dt*source;

    else
        yl = y(j - 1);
        yr = y(j + 1);
        Ql = Q(j - 1);
        Qr = Q(j + 1);

        Al = xsect.getAofY(yl);
        Ar = xsect.getAofY(yr);
        cl = xsect.getCofY(yl);
        cr = xsect.getCofY(yr);
        ul = Ql/Al;
        ur = Qr/Ar;

        Jl = [0, 1; cl^2 - ul^2, 2*ul];
        Jr = [0, 1; cr^2 - ur^2, 2*ur];

        [Fleft, spectralLeft] = Rus(yl, yj, Ql, Qj, xsect);
        [Fright, spectralRight] = Rus(yj, yr, Qj, Qr, xsect);

        AA = -0.5*lambda*(Jl + spectralLeft*I2);
        BB = I2 ...
            + 0.5*lambda*(spectralLeft + spectralRight)*I2 ...
            - dt*sourceJac;
        CC = 0.5*lambda*(Jr - spectralRight*I2);

        Aeq(row, 2*j-3:2*j+2) = [AA, BB, CC];
        Beq(row) = lambda*(Fleft - Fright) + dt*source;
    end
end
end
