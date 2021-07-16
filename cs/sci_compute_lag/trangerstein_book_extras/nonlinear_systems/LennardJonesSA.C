#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "MemoryDebugger.H"
#include "Tracer.H"
using namespace std;
int atoms=5;
int n=1;
void aFromX(const double *x,double *a) {
//Tracer t("aFromX");
  for (int i=0;i<3;i++) a[i]=0.; // atom #1 at origin
  a[3]=x[0];
  for (int i=4;i<6;i++) a[i]=0.; // atom #2 along first axis
  if (atoms>=3) {
    a[6]=x[1];
    a[7]=x[2];
    a[8]=0.; // atom #3 in 1,2 plane
  }
  if (atoms>=4) {
    for (int k=9;k<3*atoms;k++) a[k]=x[k-6];
  }
}
double lj(void *xp) {
//Tracer t("lj");
  double *x=reinterpret_cast<double*>(xp);
  double *a=new double[3*atoms];
  aFromX(x,a);
  double energy=0.;
  for (int i=0;i<atoms-1;i++) {
    int i3=i*3;
    for (int j=i+1;j<atoms;j++) {
      int j3=j*3;
      double d=pow(a[i3]-a[j3],2)+pow(a[1+i3]-a[1+j3],2)
              +pow(a[2+i3]-a[2+j3],2);
      energy+=pow(d,-6)-pow(d,-3);
    }
  }
  delete [] a; a=0;
//cout << "\tenergy = " << energy << endl;
  return energy;
}
void ljstep(const gsl_rng *r,void *xp,double step_size) {
//Tracer t("ljstep");
  double *x=reinterpret_cast<double*>(xp);
//cout << "\txp,step_size = " << xp << " " << step_size << endl;
//for (int i=0;i<n;i++) {
//  cout << "\tx[" << i << "] = " << x[i] << endl;
//}
  double norm=0.;
  for (int i=0;i<n;i++) norm+=x[i]*x[i];
  norm=sqrt(norm);
//cout << "\tnorm = " << norm << endl;
  for (int i=0;i<n;i++) {
    x[i]+=(2*gsl_rng_uniform(r)-1)*norm*step_size;
//  cout << "\tx[" << i << "] = " << x[i] << endl;
  }
}
double ljmetric(void *xp,void *yp) {
//Tracer t("ljmetric");
  double *x1=reinterpret_cast<double*>(xp);
  double *x2=reinterpret_cast<double*>(yp);
  double norm=0.;
  for (int i=0;i<n;i++) norm+=pow(x1[i]-x2[i],2);
//cout << "\tnorm = " << sqrt(norm) << endl;
  return sqrt(norm);
}
void ljprint(void *xp) {
  double *x=reinterpret_cast<double*>(xp);
  double *a=new double[3*atoms];
  aFromX(x,a);
  cout << endl;
  for (int k=0;k<atoms;k++) {
    cout << "atom[" << k << "] = <";
    for (int i=0;i<3;i++) {
      cout << a[i+k*3];
      if (i<2) cout << " , ";
    }
    cout << "> " << endl;
  }
  delete [] a; a=0;
}
void ljcopy(void *source,void *dest) {
  memcpy(dest,source,n*sizeof(double));
}
void* ljcopy_construct(void *xp) {
  double *xx=new double[n];
  memcpy(xx,xp,n*sizeof(double));
  return xx;
}
void ljdestroy(void *xp) {
  double *x=reinterpret_cast<double*>(xp);
  delete [] x; x=0;
}
void aFromX2(const int &atoms,const gsl_vector *x,double *a) {
//Tracer t("aFromX2");
  for (int i=0;i<3;i++) a[i]=0.; // atom #1 at origin
  a[3]=gsl_vector_get(x,0);
  for (int i=4;i<6;i++) a[i]=0.; // atom #2 along first axis
  if (atoms>=3) {
    a[6]=gsl_vector_get(x,1);
    a[7]=gsl_vector_get(x,2);
    a[8]=0.; // atom #3 in 1,2 plane
  }
  if (atoms>=4) {
    for (int k=9;k<3*atoms;k++) a[k]=gsl_vector_get(x,k-6);
  }
}
double lj2(const gsl_vector *x,void *params) {
//Tracer t("lj2");
  const int &atoms=*reinterpret_cast<int*>(params);
  double *a=new double[3*atoms];
  aFromX2(atoms,x,a);
  double energy=0.;
  for (int i=0;i<atoms-1;i++) {
    int i3=i*3;
    for (int j=i+1;j<atoms;j++) {
      int j3=j*3;
      double d=pow(a[i3]-a[j3],2)+pow(a[1+i3]-a[1+j3],2)
              +pow(a[2+i3]-a[2+j3],2);
      energy+=pow(d,-6)-pow(d,-3);
    }
  }
  delete [] a; a=0;
  return energy;
}
void ljgrad(const gsl_vector *x,void *params,gsl_vector *g) {
//Tracer t("ljgrad");
  const int &atoms=*reinterpret_cast<int*>(params);
  double *a=new double[3*atoms];
  aFromX2(atoms,x,a);
  double s=0.;
  double cube_root2=pow(2.,1./3.);
  double cube_root4=pow(4.,1./3.);
  for (int k=0;k<atoms;k++) {
      if (k!=1) {
      int k3=k*3;
      double d=pow(a[3]-a[k3],2)+pow(a[4]-a[1+k3],2)+pow(a[5]-a[2+k3],2);
//    double t=(pow(d,3)-2.)/pow(d,7);
      double t=(d-cube_root2)*(cube_root4+d*(cube_root2+d))/pow(d,7);
      s+=(a[3]-a[k3])*t;
    }
    gsl_vector_set(g,0,6.*s);
  }
  if (atoms>=3) {
    gsl_vector_set(g,1,0.);
    gsl_vector_set(g,2,0.);
    for (int k=0;k<atoms;k++) {
      if (k!=2) {
        int k3=k*3;
        double d=pow(a[6]-a[k3],2)+pow(a[7]-a[1+k3],2)+pow(a[8]-a[2+k3],2);
//      double t=(pow(d,3)-2.)/pow(d,7);
        double t=(d-cube_root2)*(cube_root4+d*(cube_root2+d))/pow(d,7);
        for (int i=0;i<2;i++) {
          gsl_vector_set(g,i+1,gsl_vector_get(g,i+1)+(a[i+6]-a[i+k3])*t);
        }
      }
    }
    gsl_vector_set(g,1,gsl_vector_get(g,1)*6.);
    gsl_vector_set(g,2,gsl_vector_get(g,2)*6.);
  }
  for (int j=3;j<atoms;j++) {
    int j3=j*3;
    for (int i=0;i<3;i++) gsl_vector_set(g,i+j3-6,0.);
    for (int k=0;k<atoms;k++) {
      if (k!=j) {
        int k3=k*3;
        double d=pow(a[j3]-a[k3],2)+pow(a[1+j3]-a[1+k3],2)
                +pow(a[2+j3]-a[2+k3],2);
//      double t=(pow(d,3)-2.)/pow(d,7);
        double t=(d-cube_root2)*(cube_root4+d*(cube_root2+d))/pow(d,7);
        for (int i=0;i<3;i++) {
          gsl_vector_set(g,i+j3-6,
            gsl_vector_get(g,i+j3-6)+(a[i+j3]-a[i+k3])*t);
        }
      }
    }
    for (int i=0;i<3;i++) {
      gsl_vector_set(g,i+j3-6,gsl_vector_get(g,i+j3-6)*6.);
    }
  }
  delete [] a; a=0;
}
int main(int /*argc*/,char** /*argv*/) {
  MemoryDebugger md(1);
  cout << setprecision(16);
  n=max(1,3*atoms-6);
  
  double *x=new double[n];
  x[0]=pow(2.,1./6.);
  if (atoms>=3) {
    x[1]=x[0]*0.5;
    x[2]=x[0]*sqrt(0.75);
  }
  if (atoms>=4) {
    x[3]=x[0]*0.5;
    x[4]=x[0]/sqrt(12.);
    x[5]=x[0]*sqrt(2./3.);
  }
  for (int i=6;i<n;i++) x[i]=x[i-3];
  x[n-1]=-x[n-4];
  for (int i=0;i<n;i++) {
    cout << "x[" << i << "] = " << x[i] << endl;
  }
  cout << "\tx = " << x << endl;

  double *a=new double[3*atoms];
  aFromX(x,a);

  for (int k=0;k<atoms;k++) {
    cout << "atom[" << k << "] = <";
    for (int i=0;i<3;i++) {
      cout << a[i+k*3];
      if (i<2) cout << " , ";
    }
    cout << "> " << endl;
  }
  cout << "\tenergy = " << lj(x) << endl;

  int n_tries=20;
  int iters_fixed_t=10000;
  double step_size=1.;
  double boltzman=1.;
  double initial_temperature=8.e-3;
  double damping_factor=1.003;
  double t_min=2.e-6;
  gsl_siman_params_t params={n_tries,iters_fixed_t,step_size,boltzman,
    initial_temperature,damping_factor,t_min};
  gsl_rng_env_setup();
  gsl_rng *r=gsl_rng_alloc(gsl_rng_default);

  {
    Tracer t("gsl_siman_solve");
    gsl_siman_solve(r,x,lj,ljstep,ljmetric,ljprint,ljcopy,ljcopy_construct,
      ljdestroy,sizeof(double),params);
  }

  int itmax=3000;
  gsl_vector *x2=gsl_vector_alloc(n);
  for (int i=0;i<n;i++) gsl_vector_set(x2,i,x[i]);
  gsl_multimin_function my_func2;
    my_func2.n=n;
    my_func2.f=&lj2;
    my_func2.params=reinterpret_cast<void*>(&atoms);
  gsl_vector *ss=gsl_vector_alloc(n); // step sizes
  gsl_vector_set_all(ss,1.e0);
  const gsl_multimin_fminimizer_type *T=
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s=gsl_multimin_fminimizer_alloc(T,n);
  gsl_multimin_fminimizer_set(s,&my_func2,x2,ss);
  int it=0;
  int status=GSL_CONTINUE;
  do {
    it++;
    status=gsl_multimin_fminimizer_iterate(s);
  } while (/*status==GSL_CONTINUE && */ it<itmax);
  aFromX2(atoms,s->x,a);
  for (int k=0;k<atoms;k++) {
    cout << "atom[" << k << "] = <";
    for (int i=0;i<3;i++) {
      cout << a[i+k*3];
      if (i<2) cout << " , ";
    }
    cout << "> " << endl;
  }
  cout << "\tf = " << s->fval << endl;
  gsl_vector *g2=gsl_vector_alloc(n);
  ljgrad(s->x,&atoms,g2);
  for (int i=0;i<n;i++) {
    cout << "\tgradient[" << i << "] = " << gsl_vector_get(g2,i) << endl;
  }
  gsl_vector_free(g2);
  gsl_vector_free(ss);
  gsl_vector_free(x2);
  for (int i=0;i<atoms-1;i++) {
    int i3=i*3;
    for (int j=i+1;j<atoms;j++) {
      int j3=j*3;
      cout << "distance from atom << " << i << " to atom " << j << " = "
        << sqrt(pow(a[i3]-a[j3],2)+pow(a[1+i3]-a[1+j3],2)
               +pow(a[2+i3]-a[2+j3],2)) << endl;
    }
  }

  delete [] a; a=0;
  delete [] x; x=0;
}
