clear
rate=1;
u=1;
tmax=0.9999;
nsteps=20;
tarray=zeros(1,nsteps+1);
uarray=zeros(1,nsteps+1);
dt=tmax/nsteps;
uarray(1)=u;
tarray(1)=0;
for j=1:nsteps
  uarray(j+1)=euler(rate,uarray(j),tarray(j),dt);
  tarray(j+1)=tarray(j)+dt;
end
plot(tarray,uarray)
