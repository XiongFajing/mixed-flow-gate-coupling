classdef CityGate < handle

    properties
        D
        H
        a
        Af
        Wsl
        R
        Yc
    end

    methods
        function obj = CityGate(D, H, slotRatio)
            obj.D = D;
            obj.H = H;
            obj.R = D/2;

            areaRect = D*H;
            yRect = H/2;
            areaSemi = pi*obj.R^2/2;
            ySemi = H + 4*obj.R/(3*pi);

            obj.Af = areaRect + areaSemi;
            obj.Yc = (areaRect*yRect + areaSemi*ySemi)/obj.Af;

            if nargin < 3
                slotRatio = 0.01;
            end
            obj.Wsl = obj.D*slotRatio;
            obj.a = sqrt(9.8*obj.Af/obj.Wsl);
        end

        function A = getAofY(obj, y)
            A = zeros(size(y));
            crown = obj.H + obj.R;

            for i = 1:numel(y)
                if y(i) <= 0
                    A(i) = 0;
                elseif y(i) <= obj.H
                    A(i) = obj.D*y(i);
                elseif y(i) < crown
                    hArch = y(i) - obj.H;
                    theta = 2*acos(hArch/obj.R);
                    A(i) = obj.D*obj.H ...
                        + (pi - theta + sin(theta))*obj.R^2/2;
                else
                    A(i) = obj.Af + (y(i) - crown)*obj.Wsl;
                end
            end
        end

        function w = getWofY(obj, y)
            w = zeros(size(y));
            crown = obj.H + obj.R;

            for i = 1:numel(y)
                if y(i) <= 0
                    w(i) = 0;
                elseif y(i) <= obj.H
                    w(i) = obj.D;
                elseif y(i) < crown
                    hArch = y(i) - obj.H;
                    physicalWidth = 2*sqrt(max(obj.R^2 - hArch^2, 0));

                    w(i) = max(physicalWidth, obj.Wsl);
                else
                    w(i) = obj.Wsl;
                end
            end
        end

        function P = getPofY(obj, y)
            P = zeros(size(y));
            crown = obj.H + obj.R;
            fullPerimeter = obj.D + 2*obj.H + pi*obj.R;

            for i = 1:numel(y)
                if y(i) <= 0
                    P(i) = 0;
                elseif y(i) <= obj.H
                    P(i) = obj.D + 2*y(i);
                elseif y(i) < crown
                    hArch = y(i) - obj.H;
                    theta = 2*acos(hArch/obj.R);
                    P(i) = obj.D + 2*obj.H + obj.R*(pi - theta);
                else

                    P(i) = fullPerimeter;
                end
            end
        end

        function I1 = getIofY(obj, y)
            I1 = zeros(size(y));
            crown = obj.H + obj.R;

            for i = 1:numel(y)
                if y(i) <= 0
                    I1(i) = 0;
                elseif y(i) <= obj.H
                    I1(i) = obj.D*y(i)^2/2;
                elseif y(i) < crown
                    hArch = y(i) - obj.H;
                    theta = 2*acos(hArch/obj.R);

                    areaRect = obj.D*obj.H;
                    distanceRect = y(i) - obj.H/2;
                    areaArch = (pi - theta + sin(theta))*obj.R^2/2;

                    denominator = 3*obj.R^2*(pi - theta + sin(theta));
                    centroidTerm = ( ...
                        4*obj.R^3 ...
                        - 4*max(obj.R^2 - hArch^2, 0)^(3/2))/denominator;
                    distanceArch = hArch - centroidTerm;

                    I1(i) = areaRect*distanceRect + areaArch*distanceArch;
                else
                    hSlot = y(i) - crown;
                    IFull = obj.Af*(y(i) - obj.Yc);
                    ISlot = obj.Wsl*hSlot^2/2;
                    I1(i) = IFull + ISlot;
                end
            end
        end

        function c = getCofY(obj, y)
            c = zeros(size(y));

            for i = 1:numel(y)
                if y(i) <= 0
                    c(i) = 0;
                else
                    A = obj.getAofY(y(i));
                    w = obj.getWofY(y(i));
                    c(i) = sqrt(9.8*A/max(w, obj.Wsl));
                end
            end
        end

        function y = getYofA(obj, A)
            y = zeros(size(A));
            crown = obj.H + obj.R;

            for i = 1:numel(A)
                if A(i) <= 0
                    y(i) = 0;
                elseif A(i) < obj.Af
                    lower = 0;
                    upper = crown;

                    while upper - lower > 1e-8
                        middle = 0.5*(lower + upper);
                        if obj.getAofY(middle) > A(i)
                            upper = middle;
                        else
                            lower = middle;
                        end
                    end
                    y(i) = 0.5*(lower + upper);
                else
                    y(i) = crown + (A(i) - obj.Af)/obj.Wsl;
                end
            end
        end

        function Sf = getSf(obj, nManning, y, Q)
            if isscalar(Q)
                Q = Q*ones(size(y));
            end
            Sf = zeros(size(y));

            for i = 1:numel(y)
                if y(i) <= 0
                    Sf(i) = 0;
                else
                    A = obj.getAofY(y(i));
                    P = obj.getPofY(y(i));
                    hydraulicRadius = A/P;

                    Sf(i) = nManning^2*Q(i)*abs(Q(i)) ...
                        /(A^2*hydraulicRadius^(4/3));
                end
            end
        end

        function derivative = getdSfdA(obj, nManning, y, Q)
            if isscalar(Q)
                Q = Q*ones(size(y));
            end
            derivative = zeros(size(y));

            for i = 1:numel(y)
                if y(i) <= 0
                    derivative(i) = 0;
                    continue;
                end

                A0 = obj.getAofY(y(i));
                dA = max(1e-6*max(A0, obj.Af), 1e-8);
                Aminus = max(A0 - dA, 1e-10);
                Aplus = A0 + dA;

                yMinus = obj.getYofA(Aminus);
                yPlus = obj.getYofA(Aplus);
                SfMinus = obj.getSf(nManning, yMinus, Q(i));
                SfPlus = obj.getSf(nManning, yPlus, Q(i));

                derivative(i) = (SfPlus - SfMinus)/(Aplus - Aminus);
            end
        end

        function derivative = getdSfdQ(obj, nManning, y, Q)
            if isscalar(Q)
                Q = Q*ones(size(y));
            end
            derivative = zeros(size(y));

            for i = 1:numel(y)
                if y(i) <= 0
                    derivative(i) = 0;
                else
                    A = obj.getAofY(y(i));
                    P = obj.getPofY(y(i));
                    hydraulicRadius = A/P;
                    derivative(i) = 2*nManning^2*abs(Q(i)) ...
                        /(A^2*hydraulicRadius^(4/3));
                end
            end
        end
    end
end
