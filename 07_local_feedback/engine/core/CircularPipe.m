classdef CircularPipe < handle

    properties
        D
        R
        Af
        Wsl
        a
        Yc
    end
    methods
        function obj=CircularPipe(D,slotRatio)
            validateattributes(D,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(slotRatio,{'numeric'},{'scalar','real','finite','positive'});
            obj.D=D; obj.R=D/2; obj.Af=pi*D^2/4;
            obj.Wsl=slotRatio*D; obj.a=sqrt(9.8*obj.Af/obj.Wsl); obj.Yc=D/2;
        end
        function A=getAofY(obj,y)
            A=zeros(size(y)); open=y>0 & y<obj.D; full=y>=obj.D;
            th=2*acos(max(-1,min(1,1-2*y(open)/obj.D)));
            A(open)=0.5*obj.R^2*(th-sin(th));
            A(full)=obj.Af+obj.Wsl*(y(full)-obj.D);
        end
        function w=getWofY(obj,y)
            w=zeros(size(y)); open=y>0 & y<obj.D;
            w(open)=max(2*sqrt(max(y(open).*(obj.D-y(open)),0)),obj.Wsl);
            w(y>=obj.D)=obj.Wsl;
        end
        function P=getPofY(obj,y)
            P=zeros(size(y)); open=y>0 & y<obj.D;
            P(open)=2*obj.R*acos(max(-1,min(1,1-2*y(open)/obj.D)));
            P(y>=obj.D)=pi*obj.D;
        end
        function I=getIofY(obj,y)
            I=zeros(size(y)); open=y>0 & y<obj.D; full=y>=obj.D;
            h=y(open); A=obj.getAofY(h);
            I(open)=(h-obj.R).*A+(2/3)*max(h.*(obj.D-h),0).^(3/2);
            I(full)=obj.Af*(y(full)-obj.R)+0.5*obj.Wsl*(y(full)-obj.D).^2;
        end
        function c=getCofY(obj,y)
            c=sqrt(9.8*obj.getAofY(y)./max(obj.getWofY(y),obj.Wsl));
        end
        function y=getYofA(obj,A)
            y=zeros(size(A)); full=A>=obj.Af; open=A>0 & A<obj.Af;
            y(full)=obj.D+(A(full)-obj.Af)/obj.Wsl;
            target=A(open); h=obj.D*target/obj.Af;
            lo=zeros(size(h)); hi=obj.D*ones(size(h));
            for k=1:60
                err=obj.getAofY(h)-target;
                pending=abs(err)>2e-14*obj.Af;
                if ~any(pending(:)), break; end
                above=err>0; hi(above)=h(above); lo(~above)=h(~above);
                width=2*sqrt(max(h.*(obj.D-h),0));
                candidate=h-err./max(width,eps*obj.D);
                bad=candidate<=lo | candidate>=hi | ~isfinite(candidate);
                candidate(bad)=0.5*(lo(bad)+hi(bad));
                h(pending)=candidate(pending);
            end
            y(open)=h;
        end
        function Sf=getSf(obj,n,y,Q)
            A=obj.getAofY(y); P=obj.getPofY(y);
            Sf=zeros(size(y)); wet=A>0;
            if isscalar(Q), Q=Q*ones(size(y)); end
            Sf(wet)=n^2*Q(wet).*abs(Q(wet)).*P(wet).^(4/3)./A(wet).^(10/3);
        end
        function d=getdSfdA(obj,n,y,Q)
            A=obj.getAofY(y); da=1e-6*max(A,obj.Af);
            am=max(A-da,1e-12*obj.Af); ap=A+da;
            d=(obj.getSf(n,obj.getYofA(ap),Q)-obj.getSf(n,obj.getYofA(am),Q))./(ap-am);
            d(A<=0)=0;
        end
        function d=getdSfdQ(obj,n,y,Q)
            A=obj.getAofY(y); P=obj.getPofY(y); d=zeros(size(y)); wet=A>0;
            if isscalar(Q), Q=Q*ones(size(y)); end
            d(wet)=2*n^2*abs(Q(wet)).*P(wet).^(4/3)./A(wet).^(10/3);
        end
    end
end
