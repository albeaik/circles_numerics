
clear all
nt=1e2;
nx=1e2+1;
nv=1e2;
Vmax =20; % needed when define function Vf
epsilon=1;
x=linspace(-10,10,nx);
v=linspace(0,Vmax,nv);
t=linspace(0,1,nt);

Q=zeros(nt+1,nx,nv);
Q0=@(x,v) (x>=2).*(x<=5).*(v>=10).*(v<=15);


lambdax=(t(2)-t(1))/(x(2)-x(1));
lambdav=(t(2)-t(1))/(v(2)-v(1));
[V,X]=meshgrid(v,x);
Q(1,:,:)=Q0(X,V);
Q(1,1,:)=0;
Q(1,end,:)=0;
Q(1,:,1)=0;
Q(1,:,end)=0;
d0 = 2.5; % needed when define function Vf. Change later. 
Vf=@(x) Vmax.*((tanh(x./d0-2)+tanh(2))./(1+tanh(2)));  
h=@(x) exp(-(1)./((epsilon./2).^2-(-x-epsilon/2).^2)).*(x>-epsilon).*(x<0);
theta=@(x,v) h(-x).*(Vf(x)-v).*(x<=epsilon).*(x>=0);
epsilonx=1; %viscosity parameter, needs to be changed later to a proper value
epsilonv=1;
x1 = linspace(-20,20,2.*nx);
[V,X1]=meshgrid(v,x1);
Theta=theta(X1,V);
W = zeros(nx,nv);
[V2,~]=meshgrid(v,ones(1,nx));
[V,X]=meshgrid(v,x);
for n=1:1:nt
                Q(n+1,2:nx-1,2:nv-1)=Q(n,2:nx-1,2:nv-1)-...
                lambdax.*reshape(V2(2:nx-1,2:nv-1).*reshape(((Q(n,2:nx-1,2:nv-1)+Q(n,3:nx,2:nv-1))/2-(Q(n,1:nx-2,2:nv-1)+Q(n,2:nx-1,2:nv-1))/2),nx-2,nv-2),1,nx-2,nv-2)-...
                epsilonx.*(-Q(n,3:nx,2:nv-1)+2.*Q(n,2:nx-1,2:nv-1)-Q(n,1:nx-2,2:nv-1))/2;
     n
            
            parfor i=1:1:nx
                        for j=1:1:nv
                                    W(i,j) = sum(sum(Theta(nx+i-(1:nx),j).*reshape(Q(n+1,1:nx,1:nv),nx,nv)));
                        end
            end
            W = (x(2)-x(1)).*(v(2)-v(1)).*W;
             Q(n+1,2:nx-1,2:nv-1)=Q(n+1,2:nx-1,2:nv-1)-...
             lambdav/2.*reshape(W(2:nx-1,2:nv-1).*(reshape(Q(n+1,2:nx-1,3:nv)-Q(n+1,2:nx-1,1:nv-2),nx-2,nv-2)),1,nx-2,nv-2) -...
             epsilonv.*(-Q(n+1,2:nx-1,3:nv)+2*Q(n+1,2:nx-1,2:nv-1)-Q(n+1,2:nx-1,1:nv-2))/2;
            
         figure(1)
         surf(X,V,reshape(Q(n,:,:),nx,nv));
         zlim([0,1]);
         colorbar;
         colormap(jet);
         caxis([0 1]);
         view(0,90);

         drawnow;
           % Q(n+1,2:nx-1,2:nv-1)=Q(n+1,2:nx-1,2:nv-1)-...
           % lambdav.*W(2:nx-1,2:nv-1).*(Q(n+1,2:nx-1,3:nv)-Q(n+1,2:nx-1,1:nv-2)) +...
           % epsilonv.*(-Q(n+1,2:nx-1,1:nv-2)+2*Q(n+1,2:nx-1,2:nv-1)-Q(n+1,2:nx-1,1:nv-2))/2;
end

