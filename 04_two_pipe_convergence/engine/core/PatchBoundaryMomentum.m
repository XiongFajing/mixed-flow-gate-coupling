function [M,rhs]=PatchBoundaryMomentum(M,rhs,y,qin,qout,xs,dx,dt,depth,Qold,kind)

N=numel(y); dt_dx=dt/dx;
if nargin<11,kind='direct';end
for side=1:2
    if ~isfinite(depth(side)), continue; end
    if depth(side)<=0, error('IFVM:BoundaryDepth','Prescribed depth must be positive.'); end
    if side==1, j=1; q=qin; orient=1; else, j=N; q=qout; orient=-1; end
    rows=2*j-1:2*j;
    if strcmp(kind,'riemann')
      if side==1,sgn=-1;else,sgn=1;end
      [fold,jold]=BoundaryFlux(y(j),q,xs,Qold(j),sgn,kind);
    else
      [fold,jold]=BoundaryFlux(y(j),q,xs);
    end
    [fnew,~]=BoundaryFlux(depth(side),q,xs);
    M(rows,rows)=M(rows,rows)+orient*dt_dx*jold;
    rhs(rows)=rhs(rows)+orient*dt_dx*(fnew-fold);
end
end
