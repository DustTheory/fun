#include <stdio.h> /* for printf */
#include <sys/times.h> /* for function times and struct tms */
#include <unistd.h> /* for function sysconf */
#include <values.h> /* for DBL_EPSILON */

void fcn(int *n,double *x,double *fvec,double *fjac,int *ldfjac,int *iflag) {
  if (*iflag == 1) {
    fvec[0]=10.*(x[1]-x[0]*x[0]);
    fvec[1]=1.-x[0];
  } else {
    fjac[0]=-20.*x[0];
    fjac[1]=-1.;
    fjac[*ldfjac]=10.;
    fjac[*ldfjac+1]=0.;
  }
}
      
extern void hybrj1_(void (*)(int*,double*,double*,double*,int*,int*),
  int*,double*,double*,double*,int*,double*,int*,double*,int*);

int main(int argc,char **argv) {
  int n,ldfjac,lwa,info,k;
  double elapsed,x[2],fvec[2],fjac[4],tol,wa[15];
  struct tms usage;
  clock_t start;

  n=2;
  ldfjac=2;
  lwa=15;
  tol=DBL_EPSILON;
  times(&usage);
  start=usage.tms_utime;
  for (k=0;k<4096;k++) {
    x[0]=-1.;
    x[1]=1.2;
    hybrj1_(fcn,&n,x,fvec,fjac,&ldfjac,&tol,&info,wa,&lwa);
  }
  times(&usage);
  elapsed=(double)(usage.tms_utime-start)
    /((double)(sysconf(_SC_CLK_TCK))*4096.);
  printf("average time for minpack hybrj1 = %.3g\n",elapsed);
  printf("x = %24.16g %24.16g\n",x[0],x[1]);
}
