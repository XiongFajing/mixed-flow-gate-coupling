function [M,b]=FVMMatrixFast(y,Q,n,S,dx,dt,xs,qin,qout)

y=y(:).';Q=Q(:).';N=numel(y);l=dt/dx;
A=xs.getAofY(y);c=xs.getCofY(y);u=Q./A;I=xs.getIofY(y);
Sf=xs.getSf(n,y,Q);sA=9.8*(S-Sf-A.*xs.getdSfdA(n,y,Q));
sQ=-9.8*A.*xs.getdSfdQ(n,y,Q);
sig=max(abs(u(1:end-1))+c(1:end-1),abs(u(2:end))+c(2:end));
F=[Q;Q.^2./A+9.8*I];U=[A;Q];
edge=.5*(F(:,1:end-1)+F(:,2:end)+ ...
    (U(:,1:end-1)-U(:,2:end)).*repmat(sig,2,1));
[fl,jl]=BoundaryFlux(y(1),qin,xs);[fr,jr]=BoundaryFlux(y(end),qout,xs);
rhs=l*([fl,edge]-[edge,fr]);rhs(2,:)=rhs(2,:)+dt*9.8*A.*(S-Sf);
b=rhs(:);
sl=[0,sig];sr=[sig,0];
b11=1+.5*l*(sl+sr);b12=zeros(1,N);
b21=-dt*sA;b22=b11-dt*sQ;
b12(1)=.5*l;b21(1)=b21(1)-l*jl(2,1)+.5*l*(c(1)^2-u(1)^2);
b22(1)=b22(1)+l*u(1);
b12(end)=-.5*l;b21(end)=b21(end)+l*jr(2,1)-.5*l*(c(end)^2-u(end)^2);
b22(end)=b22(end)-l*u(end);
j=1:N;rows=[2*j-1;2*j-1;2*j;2*j];cols=[2*j-1;2*j;2*j-1;2*j];
val=[b11;b12;b21;b22];
k=1:N-1;
rr=[2*k-1;2*k-1;2*k;2*k];cc=[2*k+1;2*k+2;2*k+1;2*k+2];
vv=.5*l*[-sig;ones(1,N-1);c(2:end).^2-u(2:end).^2;2*u(2:end)-sig];
rows=[rows(:);rr(:)];cols=[cols(:);cc(:)];val=[val(:);vv(:)];
rr=[2*k+1;2*k+1;2*k+2;2*k+2];cc=[2*k-1;2*k;2*k-1;2*k];
vv=-.5*l*[sig;ones(1,N-1);c(1:end-1).^2-u(1:end-1).^2;2*u(1:end-1)+sig];
M=sparse([rows;rr(:)],[cols;cc(:)],[val;vv(:)],2*N,2*N);
end
