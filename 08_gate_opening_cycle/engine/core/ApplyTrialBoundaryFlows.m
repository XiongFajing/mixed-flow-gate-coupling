function [M,b]=ApplyTrialBoundaryFlows(cache,qin,qout,dx,dt)

M=cache.M; b=cache.rhs; lambda=dt/dx; N=size(M,1);
if isfield(cache,'boundary') && strcmp(cache.boundary,'riemann')
 for side=1:2
  if side==1,ar=cache.Aleft;Q=cache.Qleft;c=cache.cleft;q=qin;j=1;ori=1;sgn=-1;
  else,ar=cache.Aright;Q=cache.Qright;c=cache.cright;q=qout;j=N/2;ori=-1;sgn=1;end
  sig0=abs(Q)/ar+c;sig=max(abs(Q),abs(2*q-Q))/ar+c;
  dF=(-2*q*Q+2*q*q)/ar+sgn*(sig*(Q-q)-sig0*Q);
  dJA=(2*q*Q-2*q*q)/ar^2;dJQ=-2*q/ar+sgn*(sig-sig0);
  M(2*j,2*j-1)=M(2*j,2*j-1)-ori*lambda*dJA;
  M(2*j,2*j)=M(2*j,2*j)-ori*lambda*dJQ;
  b(2*j-1)=b(2*j-1)+ori*lambda*q;
  b(2*j)=b(2*j)+ori*lambda*dF;
 end
 return;
end
M(2,1)=M(2,1)+lambda*(qin/cache.Aleft)^2;
M(N,N-1)=M(N,N-1)-lambda*(qout/cache.Aright)^2;
b(1:2)=b(1:2)+lambda*[qin;qin^2/cache.Aleft];
b(N-1:N)=b(N-1:N)-lambda*[qout;qout^2/cache.Aright];
end
