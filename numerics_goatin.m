clear all
close all

%Setting parameters
nt=5e3;
nx=2e2+1;
nv=2e2;
Vmax =20; % needed when define function Vf
epsilon=1;
alpha = 1;
beta = 20;
x=linspace(-0,20,nx);
v=linspace(0,20,nv);
t=linspace(0,10,nt);


%Initializing

Q=zeros(nt+1,nx,nv);
Q0=@(x,v) (x>=2).*(x<=5).*(v>=10).*(v<=15) + (x>=5).*(x<=15).*(v>=5).*(v<=10);

y = zeros(nt+1,1);
y0 = 7 ;
y(1) = y0;
w = zeros(nt+1,1);
w0 = 9*Vmax/10 ; 
w(1) = w0;

%Boundary conditions, useful functions

lambdax=(t(2)-t(1))/(x(2)-x(1));
lambdav=(t(2)-t(1))/(v(2)-v(1));
dt = t(2)-t(1);
dx = x(2) - x(1);
dv = v(2) - v(1);
[V,X]=meshgrid(v,x);
Q(1,:,:)=Q0(X,V);
Q(:,1,:)=0;
Q(:,end,:)=0;
Q(:,:,1)=0;
Q(:,:,end)=0;
d0 = 2.5; % needed when define function Vf. Change later. 
Vf=@(x) Vmax.*((tanh(x./d0-2)+tanh(2))./(1+tanh(2))); 
% h=(x>=-epsilon).*(x<=0);
h=@(x) exp(-(1)./((epsilon./2).^2-(-x-epsilon/2).^2)).*(x>-epsilon).*(x<0);

%Kernels

theta=@(x,v) alpha.*h(x).*(Vf(-x)-v).*(x>=-epsilon).*(x<=0);
theta2=@(x,v) beta.*h(x).*(-v).*(x>=-epsilon).*(x<=0);


theta3 = @(x,v,y)(alpha.*h(x-y).*(Vf(y-x)-v)); %Action of AV Bando
theta4 = @(x,v,y,w)(beta.*h(x-y).*(w-v)); %Action of AV FTL

%CFL and viscosity parameters

epsilonx=1/3; %viscosity parameter, needs to be changed later to a proper value
epsilonv=1/3;
if lambdax > 2*min(3*epsilonx,2-3*epsilonx)/7
    disp('CFL condition in x not satisfied')
elseif epsilonx>2/3 || epsilonx <0
    disp('Viscosity approximation out of range')
end

%Mesh and reshaping arrays

x1 = linspace(-20,20,2.*nx);
[V,X1] = meshgrid(v,x1);
Theta = theta(X1,V);
Thetaflip = flip(Theta);
%We want that Theta(i-l-1,j)=Thetab(i,l,j)
%see below the previous line for W that is now commented, and the definition of Thetab below to check that we
%have the same thing.

Thetab = zeros(nx,nv,nx);
for i = 1:nx
    Thetab(i,1:nv,1:nx) = (Thetaflip(nx-i+(1:nx),1:nv))';
end

%COMMENT NEXT IF DON'T WANT H2
v1 = linspace(0,20,2.*nv);
% [V1,X11] = meshgrid(v1,x1);
% Theta2 = theta2(X11,V1);
% Thetaflip20 = flip(Theta2);
% Thetaflip2 = flip(Thetaflip20,2);

%LIMITING FACTOR
% Theta2b = zeros(nx,nv,nx,nv);
% Theta2b0 = zeros(nx,nx,2*nv);
% for i = 1:nx
%     for j = 1:nv
%         Theta2b(i,j,1:nv,1:nx) = (Thetaflip2(nx-i+(1:nx),nv-j+(1:nv)))';
%     end
% end
% 
% for i = 1:nx
% %    for j = 1:nv
%         Theta2b0(i,1:nx,1:2*nv) = Thetaflip2(nx-i+(1:nx),(1:2*nv));
% %    end
% end
% 
% for j = 1:nv
% %    for j = 1:nv
%         Theta2b(1:nx,j,1:nx,1:nv) = Theta2b0(1:nx,1:nx,nv-j+(1:nv));
%         %(Thetaflip2(nx-i+(1:nx),(1:2*nv)));
% %    end
% end

%------ H2--------
Thetax = h(fliplr(x1));
Thetav = fliplr(v1);
Thetaxb = zeros(nx,nx);
Thetavb = zeros(nv,nv);
for i = 1:nx
    Thetaxb(i,1:nx) = (Thetax(nx-i+(1:nx)));
end

for j = 1:nv
    Thetavb(j,1:nv) = (Thetav(nv-j+(1:nv)));
end

W = zeros(nx,nv);
[V2,~]=meshgrid(v,ones(1,nx));
[V,X]=meshgrid(v,x);

%Time iteration


for n=1:1:nt
                Q(n+1,2:nx-1,2:nv-1)=Q(n,2:nx-1,2:nv-1)-...
                lambdax.*reshape(V2(2:nx-1,2:nv-1).*reshape(((Q(n,2:nx-1,2:nv-1)+Q(n,3:nx,2:nv-1))/2-(Q(n,1:nx-2,2:nv-1)+Q(n,2:nx-1,2:nv-1))/2),nx-2,nv-2),1,nx-2,nv-2)-...
                epsilonx.*(-Q(n,3:nx,2:nv-1)+2.*Q(n,2:nx-1,2:nv-1)-Q(n,1:nx-2,2:nv-1))/2;
     n
            %Q(n+1,2:nx-1,2:nv-1)=Q(n,2:nx-1,2:nv-1)-...
            %lambdax.*v(2:nv-1).*((Q(n,2:nx-1,2:nv-1)+Q(n,3:nx,2:nv-1))/2-(Q(n,1:nx-2,2:nv-1)+Q(n,2:nx-1,2:nv-1))/2)+...
            %epsilonx.*(-Q(n,3:nx,2:nv-1)+2.*Q(n,2:nx-1,2:nv-1)-Q(n,1:nx-2,2:nv-1))/2;
            %n
            
  %          for i=1:1:nx
                      % for j=1:1:nv
                       %             W(i,j) = sum(sum(Theta(nx+i-(1:nx),j).*reshape(Q(n+1,1:nx,1:nv),nx,nv)));
                       % end
                       %donc Thetaflip(nx-i+(1:nx),:)=Theta(nx-i+(1:nx),:)
                                    
                       %W(i,1:nv) = sum(Thetaflip(nx-i+(1:nx),:)'*reshape(Q(n+1,1:nx,1:nv),nx,nv),2);%            W(1:nx,1:nv) = sum(einsum(Thetab(1:nx,1:nv,1:nx),reshape(Q(n+1,1:nx,1:nv),nx,nv),3,1),3);
            
%COMMENT NEXT LINE AND UNCOMMENT PREVIOUS ONE IF DON'T WANT H2 - NOT CHECKED, HANDLE CAREFULLY
%W(1:nx,1:nv) = sum(einsum(Thetab(1:nx,1:nv,1:nx),reshape(Q(n+1,1:nx,1:nv),nx,nv),3,1),3) + sum(einsum(Theta2b,reshape(Q(n+1,1:nx,1:nv),nx,nv),[3 4],[1 2]))+ theta3(X,V,y(n)) + theta4(X,V,y(n),w(n));
W(1:nx,1:nv) = sum(einsum(Thetab(1:nx,1:nv,1:nx),reshape(Q(n+1,1:nx,1:nv),nx,nv),3,1),3) - (beta * Thetavb * (Thetaxb * reshape(Q(n+1,1:nx,1:nv),nx,nv))')' + theta3(X,V,y(n)) + theta4(X,V,y(n),w(n));
%sum(dot(Thetaflip2(nx-:+(1:nx),nv-:+(1:nv)),reshape(Q(n+1,1:nx,1:nv),nx,nv)))
%            end
            W = (x(2)-x(1)).*(v(2)-v(1)).*W;
             Q(n+1,2:nx-1,2:nv-1)=Q(n+1,2:nx-1,2:nv-1)-...
             lambdav/2.*reshape(W(2:nx-1,2:nv-1).*(reshape(Q(n+1,2:nx-1,3:nv)-Q(n+1,2:nx-1,1:nv-2),nx-2,nv-2)),1,nx-2,nv-2) -...
             epsilonv.*(-Q(n+1,2:nx-1,3:nv)+2*Q(n+1,2:nx-1,2:nv-1)-Q(n+1,2:nx-1,1:nv-2))/2;
            %AV ODE update
            y(n+1) = y(n) + w(n) * dt;
            w(n+1) = w(n) + W(sum((y(n)-x)>=0),sum((w(n)-v)>=0)) * dt;
            y(n+1) = round(y(n+1)./dx).*dx;
            w(n+1) = round(w(n+1)./dv).*dv;

         figure(1)
         
         surf(X,V,reshape(Q(n,:,:),nx,nv));
         shading interp
         title(['time= ',num2str(t(n))])
         zlim([0,1]);
         colorbar;
         colormap(jet);
         caxis([0 1]);
         view(0,90);
         hold on
         plot(y(n+1),w(n+1),'r*')
         hold off

         drawnow;
           % Q(n+1,2:nx-1,2:nv-1)=Q(n+1,2:nx-1,2:nv-1)-...
           % lambdav.*W(2:nx-1,2:nv-1).*(Q(n+1,2:nx-1,3:nv)-Q(n+1,2:nx-1,1:nv-2)) +...
           % epsilonv.*(-Q(n+1,2:nx-1,1:nv-2)+2*Q(n+1,2:nx-1,2:nv-1)-Q(n+1,2:nx-1,1:nv-2))/2;
end
