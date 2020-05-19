nt=1e2;
nx=1e3;
nv=1e3;
Vmax =20; % needed when define function Vf
x=linspace(-10,10,nx);
v=linspace(0,Vmax,nv);
t=linspace(0,10,nt);

Q=zeros(nt+1,nx,nv);
lambdax=(t(2)-t(1))/(x(2)-x(1));
lambdav=(t(2)-t(1))/(v(2)-v(1));
Q(1,:,:)=1/2*ones(nx,nv);
d0 = 2.5; % needed when define function Vf. Change later. 
Vf=@(x) Vmax.*((tanh(x./d0-2)+tanh(2))./(1+tanh(2)));  
h=@(x) exp(-(1)./((epsilonx./2).^2-(-x-epsilonx/2).^2)).*(x>-epsilonx).*(x<0);
theta=@(x,v) h(-x).*(Vf(x)-v).*(x<=epsilonx).*(x>=0);
epsilonx=1; %viscosity parameter, needs to be changed later to a proper value
epsilonv=1;
x1 = linspace(-20,20,2.*nx);
[V,X1]=meshgrid(v,x1);
Theta=theta(X1,V);
W = zeros(nx,nv);

for n=1:1:nt
                Q(n+1,2:nx-1,2:nv-1)=Q(n,2:nx-1,2:nv-1)-...
                lambdax.*v(2:nv-1).*((Q(n,2:nx-1,2:nv-1)+Q(n,3:nx,2:nv-1))/2-(Q(n,1:nx-2,2:nv-1)+Q(n,2:nx-1,2:nv-1))/2)+...
                epsilonx.*(-Q(n,3:nx,2:nv-1)+2.*Q(n,2:nx-1,2:nv-1)-Q(n,1:nx-2,2:nv-1))/2;
     
            %Q(n+1,2:nx-1,2:nv-1)=Q(n,2:nx-1,2:nv-1)-...
            %lambdax.*v(2:nv-1).*((Q(n,2:nx-1,2:nv-1)+Q(n,3:nx,2:nv-1))/2-(Q(n,1:nx-2,2:nv-1)+Q(n,2:nx-1,2:nv-1))/2)+...
            %epsilonx.*(-Q(n,3:nx,2:nv-1)+2.*Q(n,2:nx-1,2:nv-1)-Q(n,1:nx-2,2:nv-1))/2;
            %n
            
            for i=1:1:nx
                        for j=1:1:nv
                                    s = zeros(1,nv);
                                    for l = -nx:nx
%                                                s=s+Theta(nx+i-l,j).*Q(n+1,l+nx,2:nv-1);
                                                 s=s+Theta(nx+i-l,j).*Q(n+1,l,1:nv); %Theta(1) corresponds in fact to -nx
                                   end
                                    W(i,j) = sum(s);
                        end
            end
            W = (x(2)-x(1)).*(v(2)-v(1)).*W;
             Q(n+1,2:nx-1,2:nv-1)=Q(n+1,2:nx-1,2:nv-1)-...
             lambdav.*W(2:nx-1,2:nv-1).*(Q(n+1,2:nx-1,3:nv)-Q(n+1,2:nx-1,1:nv-2)) +...
             epsilonv.*(-Q(n+1,2:nx-1,3:nv-1)+2*Q(n+1,2:nx-1,2:nv-1)-Q(n+1,2:nx-1,1:nv-2))/2;
            
           % Q(n+1,2:nx-1,2:nv-1)=Q(n+1,2:nx-1,2:nv-1)-...
           % lambdav.*W(2:nx-1,2:nv-1).*(Q(n+1,2:nx-1,3:nv)-Q(n+1,2:nx-1,1:nv-2)) +...
           % epsilonv.*(-Q(n+1,2:nx-1,1:nv-2)+2*Q(n+1,2:nx-1,2:nv-1)-Q(n+1,2:nx-1,1:nv-2))/2;
end
