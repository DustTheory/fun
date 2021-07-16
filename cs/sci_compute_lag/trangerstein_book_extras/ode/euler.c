#include <stdio.h> /* for FILE, scanf, printf, fprintf */
#include <stdlib.h> /* for EXIT_SUCCESS */

/*main program must have the following argument list and return value*/
int main(int argc,char **argv) {
/*must declare all variables first*/
  int i=0,nsteps=0; 
  float dt=1.,rate=1.,t=0.,tmax=1.,u=1.;
  FILE *out;

/*function definitions*/
  float f(float u,float t) { return rate*u; }
  float euler(float u,float t,float dt) { return u+dt*f(u,t); }

/*read input from stdin*/
  printf("enter rate\n");
  scanf("%f",&rate);
  printf("enter initial u\n");
  scanf("%f",&u);
  printf("enter tmax\n");
  scanf("%f",&tmax);
  printf("enter number steps\n");
  scanf("%d",&nsteps);
  printf("rate = %f u = %f\n",rate,u);
  printf("tmax = %f number steps = %d\n",tmax,nsteps);

/*apply forward-Euler method*/
  t=0.;
  dt=tmax/((float) nsteps);
  out=fopen("coutput","w");
  fprintf(out,"%f %f\n",t,u);
  for (i=1;i<=nsteps;i++) {
    u=euler(u,t,dt);
    t+=dt;
    fprintf(out,"%f %f\n",t,u);
  }
  fclose(out);
  return EXIT_SUCCESS;
}
