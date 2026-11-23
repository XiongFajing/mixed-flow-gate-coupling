function G=TrajOpening(t,G0,e,tc,tr,to)
G=G0(:);
if t<=0,G(end)=.1;
elseif t<tr
 if tc>0,G(end)=.1*max(1-t/tc,0);else,G(end)=0;end
else
 if to>0,G(end)=e*min((t-tr)/to,1);else,G(end)=e;end
end
end
