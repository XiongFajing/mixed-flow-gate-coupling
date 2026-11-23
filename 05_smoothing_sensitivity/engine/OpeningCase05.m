function G=OpeningCase05(t,G0,u)
f=ones(3,1);
f(1)=1-min(max((t-u.endGateStart)/u.endGateDuration(1),0),1);
f(3)=1-min(max((t-u.endGateStart)/u.endGateDuration(2),0),1);

if t<=u.middleTimes(1)
    f(2)=u.middleFactors(1);
elseif t>=u.middleTimes(end)
    f(2)=u.middleFactors(end);
else
    f(2)=interp1(u.middleTimes,u.middleFactors,t,'linear');
end
G=G0(:).*f;
end
