function G=OpeningCase05(t,G0,u)
f=ones(3,1);
f(1)=1-min(max((t-u.endGateStart)/u.endGateDuration(1),0),1);
f(3)=1-min(max((t-u.endGateStart)/u.endGateDuration(2),0),1);
f(2)=1-(1-u.middleGateFinalFactor)*min(max((t-u.middleGateStart)/u.middleGateDuration,0),1);
G=G0(:).*f;
end
