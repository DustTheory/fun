function w=euler(rate,u,t,dt)
  w=u+dt*fode(rate,u,t);
