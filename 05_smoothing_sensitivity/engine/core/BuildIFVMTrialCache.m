function cache=BuildIFVMTrialCache(y,Q,nManning,bedSlope,dx,dt,xsect,numerics)

if nargin<8,numerics.kernel='rusanov';end
if ~isfield(numerics,'boundary'),numerics.boundary='direct';end
if ~strcmp(numerics.kernel,'rusanov')
 error('IFVM:Kernel','This release preserves the uploaded Rusanov interior kernel.');
end
if ~any(strcmp(numerics.boundary,{'direct','riemann'})),error('Unknown boundary momentum closure.');end
cache=cell(numel(y),1);
for p=1:numel(y)
    [M,b]=FVMMatrixFast(y{p},Q{p},nManning(p),bedSlope(p),dx,dt,xsect(p),0,0);
    if strcmp(numerics.boundary,'riemann')
      for side=1:2
        if side==1,j=1;ori=1;sgn=-1;else,j=numel(y{p});ori=-1;sgn=1;end
        rows=2*j-1:2*j;
        [fold,jold]=BoundaryFlux(y{p}(j),0,xsect(p));
        [fnew,jnew]=BoundaryFlux(y{p}(j),0,xsect(p),Q{p}(j),sgn,'riemann');
        M(rows,rows)=M(rows,rows)+ori*dt/dx*(jold-jnew);
        b(rows)=b(rows)+ori*dt/dx*(fnew-fold);
      end
    end
    cache{p}=struct('M',M,'rhs',b,'Aleft',xsect(p).getAofY(y{p}(1)), ...
        'Aright',xsect(p).getAofY(y{p}(end)));
    cache{p}.boundary=numerics.boundary;
    cache{p}.Qleft=Q{p}(1);cache{p}.Qright=Q{p}(end);
    cache{p}.cleft=xsect(p).getCofY(y{p}(1));cache{p}.cright=xsect(p).getCofY(y{p}(end));
end
end
