clear all
close all

%Setting parameters
nt=3e4;
nx=0.5e2+1;
nv=0.5e2;
Vmax =20; % needed when define function Vf
epsilon=2;
alpha = 1;
beta = 20;
x=linspace(-0,20,nx);
v=linspace(0,20,nv);
t=linspace(0,10,nt);


%Initializing

Q=zeros(nt+1,nx,nv);
Q0=@(x,v) (x>=2).*(x<=5).*(v>=10).*(v<=15) + (x>=5).*(x<=15).*(v>=5).*(v<=10);

M = 2; %number of AVs
y = zeros(nt+1,M);
y0 = [7;9];
y(1,:) = y0;
w = zeros(nt+1,M);
w0 = [9*Vmax/10;9*Vmax/10]; 
w(1,:) = w0;

%Boundary conditions, useful functions

lambdax=(t(2)-t(1))/(x(2)-x(1));
lambdav=(t(2)-t(1))/(v(2)-v(1));
dt = t(2)-t(1);
dx = x(2) - x(1);
dv = v(2) - v(1);
%[V,X]=meshgrid(v,x);
Q(1,:,:)=rand(nx,nv);
Q(:,1,:)=0;
Q(:,end,:)=0;
Q(:,:,1)=0;
Q(:,:,end)=0;
d0 = 2.5; % needed when define function Vf. Change later. 
Vf=@(x) Vmax.*((tanh(x./d0-2)+tanh(2))./(1+tanh(2))); 
% h=@(x)1/epsilon.*(x>=-epsilon).*(x<=0);
%h=@(x) exp(-(1)./((epsilon./2).^2-(-x-epsilon/2).^2)).*(x>-epsilon).*(x<0);
%h=@(x) max(exp(-(1)./((epsilon./2).^2-(-x-epsilon/2).^2)).*(x>=-epsilon).*(x<=0), 0,'omitnan');
%epsilon_0 = 1
%epsilon_1 = 2
c = epsilon_0.^2.*(epsilon_0+epsilon_1)/(32.*exp(2));
h0 = @(x) c.^(-1).*max((x.^2).*exp(-x.*4./epsilon_0).*(x<=epsilon_0./2).*(x>=0)+...
(epsilon_0.^2)./(epsilon_1.^2).*(((epsilon_1+epsilon_0)./2-x).^2).*exp(-((epsilon_1+epsilon_0)./2-x).*4./epsilon_1).*(x<=(epsilon_1+epsilon_0)./2).*(x>epsilon_0./2), 0,'omitnan');
h = @(x)(h0(-x));
plot(-(0:0.01:2),h(-(0:0.01:2)))

%Kernels

%theta=@(x,v) alpha.*h(x).*(Vf(-x)-v).*(x>=-epsilon).*(x<=0);
%theta2=@(x,v) beta.*h(x).*(-v).*(x>=-epsilon).*(x<=0);

theta1x=@(x)(alpha.*(h(x)./x.^2));

theta1av=@(v)((T+v./(2*sqrt(ab))).^2);
theta2av=@(v)(2.*s0.*(T+v./(2*sqrt(ab))));
theta3av=@(v)(s_0.^2);

theta1a=@(x,v)(theta1x(x).*theta1av(v));
theta2a=@(x,v)(theta1x(x).*theta2av(v));
theta3a=@(x,v)(theta1x(x).*theta3av(v));

%theta1a=@(x,v)(theta1x(x).*theta1av(v));alpha.*(h(x)./x.^2).*(T+v./(2*sqrt(ab))).^2);
%theta2a=@(x,v)(alpha.*2.*(h(x)./x.^2).*s0.*(T+v./(2*sqrt(ab)))); 
%theta3a=@(x,v)(alpha.*(h(x)./x.^2).*s0.^2);


theta3 = @(x,v,y)(alpha.*h(x-y).*(Vf(y-x)-v)./M); %Action of AV Bando
theta4 = @(x,v,y,w)(beta.*h(x-y).*(w-v)./M); %Action of AV FTL

%CFL and viscosity parameters

epsilonx=1/3; %viscosity parameter, needs to be changed later to a proper value
epsilonv=1/3;
if lambdax > 2*min(3*epsilonx,2-3*epsilonx)/(1+6*v(end))
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

%-------New H, IDM--------
Thetax = theta1x(fliplr(x1));
Theta1v = theta1av(fliplr(v1));
Theta2v = theta2av(fliplr(v1));
Theta3v = theta3av(fliplr(v1));

Thetaxb = zeros(nx,nx);
Thetavb = zeros(nv,nv);
for i = 1:nx
%    Thetaxb(i,1:nx) = (Thetax(nx-i+(1:nx)));
end

for j = 1:nv
%    Thetavb(j,1:nv) = (Thetav(nv-j+(1:nv)));
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
ThetaW3M = zeros(nx, nv);
ThetaW4M = zeros(nx, nv);
for k = 1:M
    ThetaW3M =  ThetaW3M + theta3(X,V,y(n,k));
    ThetaW4M =  ThetaW4M + theta4(X,V,y(n,k),w(n,k));
end
W(1:nx,1:nv) = sum(einsum(Thetab(1:nx,1:nv,1:nx),reshape(Q(n+1,1:nx,1:nv),nx,nv),3,1),3) - (beta * Thetavb * (Thetaxb * reshape(Q(n+1,1:nx,1:nv),nx,nv))')'+ ThetaW3M + ThetaW4M;
%sum(dot(Thetaflip2(nx-:+(1:nx),nv-:+(1:nv)),reshape(Q(n+1,1:nx,1:nv),nx,nv)))
%            end
            W = (x(2)-x(1)).*(v(2)-v(1)).*W;
             Q(n+1,2:nx-1,2:nv-1)=Q(n+1,2:nx-1,2:nv-1)-...
             lambdav/2.*reshape(W(2:nx-1,2:nv-1).*(reshape(Q(n+1,2:nx-1,2:nv-1)+Q(n+1,2:nx-1,3:nv),nx-2,nv-2)),1,nx-2,nv-2) +...
             lambdav/2.*reshape(W(2:nx-1,1:nv-2).*(reshape(Q(n+1,2:nx-1,2:nv-1)+Q(n+1,2:nx-1,1:nv-2),nx-2,nv-2)),1,nx-2,nv-2) -...
             epsilonv.*(-Q(n+1,2:nx-1,3:nv)+2*Q(n+1,2:nx-1,2:nv-1)-Q(n+1,2:nx-1,1:nv-2))/2;
         if lambdav > 2*min(3*epsilonv,2-3*epsilonv)/(1+6*max(max(Q(n+1,:,:))))
                    disp('CFL condition in y not satisfied')
         elseif epsilonv>2/3 || epsilonv <0
                    disp('Viscosity approximation out of range')
        end
         
         
            %AV ODE update
            for k = 1:M
                y(n+1,k) = y(n,k) + w(n,k) * dt;
                w(n+1,k) = w(n,k) + W(sum((y(n,k)-x)>=0),sum((w(n,k)-v)>=0)) * dt;
            end
            %y(n+1) = round(y(n+1)./dx).*dx;
            %w(n+1) = round(w(n+1)./dv).*dv;

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
         for k = 1:M
            plot3(y(n+1,k),w(n+1,k),1,'r*')
         end
         hold off

         drawnow;
           % Q(n+1,2:nx-1,2:nv-1)=Q(n+1,2:nx-1,2:nv-1)-...
           % lambdav.*W(2:nx-1,2:nv-1).*(Q(n+1,2:nx-1,3:nv)-Q(n+1,2:nx-1,1:nv-2)) +...
           % epsilonv.*(-Q(n+1,2:nx-1,1:nv-2)+2*Q(n+1,2:nx-1,2:nv-1)-Q(n+1,2:nx-1,1:nv-2))/2;
end
