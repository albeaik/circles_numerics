nt=10000;
nx=10000;
dt=0.01;
dx=0.01;

t=dt.*(1:1:nt);
x=dx.*(-nx-0.5:1:nx+0.5);
v=dv.*(-nx-0.5:1:nx+0.5);

lx=dt./dx;
lv=dt./dv;
hallo
