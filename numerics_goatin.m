nt=1e2;
nx=1e3;
nv=1e3;
x=linspace(-10,10,nx);
v=linspace(0,20,nv);
t=linspace(0,10,nt);

Q=zeros(nt+1,nx,nv);
lambdax=(t(2)-t(1))/(x(2)-x(1));
lambdav=(t(2)-t(1))/(v(2)-v(1));
Q(1,:,:)=1/2*ones(nx,nv);
V=@(x) x;   %needs to be changed to the actual functions
h=@(x) x;
theta=@(x,v) h(-x).*(V(x)-v).*(x<=epsilon).*(x>=0);
epsilonx=1; %viscosity parameter, needs to be changed later to a proper value
epsilonv=1;
[X,V]=meshgrid(x,v);
Theta=theta(X,V);

for n=1:1:nt
            Q(n+1,2:nx-1,2:nv-1)=Q(n,2:nx-1,2:nv-1)-lambdax.*v(2:nv-1).*((Q(n,2:nx-1,2:nv-1)+Q(n,3:nx,2:nv-1))/2-(Q(n,1:nx-2,2:nv-1)+Q(n,2:nx-1,2:nv-1))/2)+epsilonx.*(-Q(n,3:nx,2:nv-1)+2.*Q(n,2:nx-1,2:nv-1)-Q(n,1:nx-2,2:nv-1))/2;
            n
            
            for i=1:1:nx
                        for j=1:1:nv
                        W(i,j)=(x(2)-x(1)).*(v(2)-v(1))*(Theta(i,
                        end
            end
            
            W=?????
            Q(n+1,2:nx-1,2:nv-1)=Q(n+1,2:nx-1,2:nv-1)-lambdav.*W(2:nx-1,2:nv-1).*(Q(n+1,2:nx-1,3:nv)-Q(n+1,2:nx-1,3:nv)) +epsilonv.*(-Q(n+1,2:nx-1,1:nv-2)+2*Q(n+1,2:nx-1,2:nv-1)-Q(n+1,2:nx-1,1:nv-2))/2;
end
