function q = GateFlow(Hu,Hd,opening,Cd,xsect,headEps,allowReverse,gateLaw)

if strcmp(gateLaw,'smooth')
    q=SmoothGateFlow(Hu,Hd,opening,Cd,xsect,headEps,allowReverse);
elseif strcmp(gateLaw,'classical')
    if opening<=0, q=0; return; end
    area=min(xsect.getAofY(opening),xsect.Af);
    dH=Hu-Hd;
    if allowReverse
        q=Cd*area*sqrt(2*9.8)*sign(dH)*sqrt(abs(dH));
    else
        q=Cd*area*sqrt(2*9.8)*sqrt(max(dH,0));
    end
else
    error('IFVM:GateLaw','gateLaw must be classical or smooth.');
end
end
