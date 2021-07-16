#include <cstdlib> // for EXIT_SUCCESS
#include <iostream> // for cout
#include <limits> // for numeric_limits
#include <sys/times.h> // for function times and struct tms
#include <unistd.h> // for function sysconf

using namespace std;

void fcn(const int &n,const double *x,double *fvec,double *fjac,
const int &ldfjac,const int &iflag) {
  if (iflag == 1) {
    fvec[0]=10.*(x[1]-x[0]*x[0]);
    fvec[1]=1.-x[0];
  } else {
    fjac[0]=-20.*x[0];
    fjac[1]=-1.;
    fjac[ldfjac]=10.;
    fjac[ldfjac+1]=0.;
  }
}
      
extern "C" {
  void hybrj1_(
    void (*)(const int&,const double*,double*,double*,const int&,const int&),
    const int&,double*,double*,double*,const int&,const double&,int&,
    double*,const int&);
}

int main(int /*argc*/,char** /*argv*/) {
  int info;
  double x[2],fvec[2],fjac[4],tol,wa[15];
  struct tms usage;
  clock_t start;

  tol=numeric_limits<double>::epsilon();
  times(&usage);
  start=usage.tms_utime;
  for (int k=0;k<4096;k++) {
    x[0]=-1.;
    x[1]=1.2;
    hybrj1_(fcn,2,x,fvec,fjac,2,tol,info,wa,15);
  }
  times(&usage);
  double elapsed=(double)(usage.tms_utime-start)
    /((double)(sysconf(_SC_CLK_TCK))*4096.);
  cout << "average time for minpack hybrj1 = " << elapsed << endl;
  cout << "x = " << x[0] << " " << x[1] << endl;
  return EXIT_SUCCESS;
}
