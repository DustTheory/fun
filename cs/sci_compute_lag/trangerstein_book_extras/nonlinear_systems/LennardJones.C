#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "MemoryDebugger.H"
#include "Tracer.H"
using namespace std;
void aFromX(const int &atoms,const gsl_vector *x,double *a) {
//Tracer t("aFromX");
//cout << "\tatoms,x,a = " << atoms << " " << x << " " << a << endl;
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
//for (int k=0;k<atoms;k++) {
//  cout << "atom[" << k << "] = <";
//  int k3=k*3;
//  for (int i=0;i<3;i++) {
//    cout << a[i+k3];
//    if (i<2) cout << " , ";
//  }
//  cout << "> " << endl;
//}
}
double lj(const gsl_vector *x,void *params) {
//Tracer t("lj");
  const int &atoms=*reinterpret_cast<int*>(params);
  double *a=new double[3*atoms];
  aFromX(atoms,x,a);
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
  aFromX(atoms,x,a);
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
void ljfdf(const gsl_vector *x,void *params,double *f,gsl_vector *g) {
//Tracer t("ljfdf");
  const int &atoms=*reinterpret_cast<int*>(params);
  double *a=new double[3*atoms];
  aFromX(atoms,x,a);
  double energy=0.;
  double cube_root2=pow(2.,1./3.);
  double cube_root4=pow(4.,1./3.);

//i=0:
  for (int j=1;j<atoms;j++) {
    int j3=j*3;
    double d=pow(a[0]-a[j3],2)+pow(a[1]-a[1+j3],2) +pow(a[2]-a[2+j3],2);
    energy+=pow(d,-6)-pow(d,-3);
  }

//i=1;
  double s=0.;
  for (int j=0;j<atoms;j++) {
    if (j!=1) {
      int j3=j*3;
      double d=pow(a[3]-a[j3],2)+pow(a[4]-a[1+j3],2)+pow(a[5]-a[2+j3],2);
      energy+=pow(d,-6)-pow(d,-3);
//    double t=(pow(d,3)-2.)/pow(d,7);
      double t=(d-cube_root2)*(cube_root4+d*(cube_root2+d))/pow(d,7);
      s+=(a[3]-a[j3])*t;
    }
    gsl_vector_set(g,0,6.*s);
  }
  if (atoms>=3) {
//  i=2;
    gsl_vector_set(g,1,0.);
    gsl_vector_set(g,2,0.);
    for (int j=0;j<atoms;j++) {
      if (j!=2) {
        int j3=j*3;
        double d=pow(a[6]-a[j3],2)+pow(a[7]-a[1+j3],2)+pow(a[8]-a[2+j3],2);
        energy+=pow(d,-6)-pow(d,-3);
//      double t=(pow(d,3)-2.)/pow(d,7);
        double t=(d-cube_root2)*(cube_root4+d*(cube_root2+d))/pow(d,7);
        for (int i=0;i<2;i++) {
          gsl_vector_set(g,i+1,gsl_vector_get(g,i+1)+(a[i+6]-a[i+j3])*t);
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
        energy+=pow(d,-6)-pow(d,-3);
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
  *f=energy*0.5;
  delete [] a; a=0;
}
int main(int /*argc*/,char** /*argv*/) {
  MemoryDebugger md(1);
  cout << setprecision(16);
  int atoms=4;
  int n=max(1,3*atoms-6);
  int itmax=3000;
//cout << "atoms,n,itmax = " << atoms << " " << n << " " << itmax << endl;
  
  gsl_vector *x=gsl_vector_alloc(n);
  gsl_vector *xx=gsl_vector_alloc(n);
  gsl_vector_set(x,0,pow(2.,1./6.));
//gsl_vector_set(x,0,1.);
  if (atoms>=3) {
    gsl_vector_set(x,1,gsl_vector_get(x,0)*0.5);
    gsl_vector_set(x,2,gsl_vector_get(x,0)*sqrt(0.75));
  }
  if (atoms>=4) {
    gsl_vector_set(x,3,gsl_vector_get(x,0)*0.5);
    gsl_vector_set(x,4,gsl_vector_get(x,0)/sqrt(12.));
    gsl_vector_set(x,5,gsl_vector_get(x,0)*sqrt(2./3.));
  }
  for (int i=0;i<n;i++) {
    cout << "x[" << i << "] = " << gsl_vector_get(x,i) << endl;
  }

  double *a=new double[3*atoms];
  aFromX(atoms,x,a);

  for (int k=0;k<atoms;k++) {
    cout << "atom[" << k << "] = <";
    for (int i=0;i<3;i++) {
      cout << a[i+k*3];
      if (i<2) cout << " , ";
    }
    cout << "> " << endl;
  }
  cout << "\tenergy = " << lj(x,&atoms) << endl;
  gsl_vector *ggg=gsl_vector_alloc(n);
  ljgrad(x,&atoms,ggg);
  for (int i=0;i<n;i++) {
    cout << "\tgradient[" << i << "] = "
         << gsl_vector_get(ggg,i) << endl;
  }
  gsl_vector_free(ggg);

  double step_size=0.1;
  double tol=pow(numeric_limits<double>::epsilon(),1./3.);

/*
  gsl_multimin_function_fdf my_func;
  my_func.n=n;
  my_func.f=&lj;
  my_func.df=&ljgrad;
  my_func.fdf=&ljfdf;
  my_func.params=reinterpret_cast<void*>(&atoms);
  for (int it2=0;it2<20;it2++) {
    for (int i=0;i<n;i++) {
      gsl_vector_set(xx,i,gsl_vector_get(x,i)*(1.+drand48()*1.e1));
//    gsl_vector_set(xx,i,gsl_vector_get(x,i));
    }
    const gsl_multimin_fdfminimizer_type *T=
      gsl_multimin_fdfminimizer_vector_bfgs2;
    gsl_multimin_fdfminimizer *s=gsl_multimin_fdfminimizer_alloc(T,n);
    gsl_multimin_fdfminimizer_set(s,&my_func,xx,step_size,tol);
    int it=0;
//  see gsl_errno.h for status codes
    int status=GSL_CONTINUE;
    do {
//    Tracer t("iteration loop");
      it++;
      status=gsl_multimin_fdfminimizer_iterate(s);
//    cout << "\tminimizer status = " << status << endl;
//    if (status!=GSL_CONTINUE) break;
      status=gsl_multimin_test_gradient(s->gradient,abs(s->f)*1.e-3);
//    cout << "\ttest gradients status = " << status << endl;
//    for (int i=0;i<n;i++) {
//      cout << "\tx[" << i << "] = " << gsl_vector_get(s->x,i) << endl;
//    }
//    cout << "\tf = " << s->f << endl;
//    for (int i=0;i<n;i++) {
//      cout << "\tgradient[" << i << "] = "
//           << gsl_vector_get(s->gradient,i) << endl;
//    }
    } while (status==GSL_CONTINUE && it<itmax);
    cout << "\tstatus = " << status << endl;
//  if (status==GSL_SUCCESS) {
      cout << "\ttrial = " << it2 << endl;
      aFromX(atoms,s->x,a);
      for (int k=0;k<atoms;k++) {
        cout << "atom[" << k << "] = <";
        for (int i=0;i<3;i++) {
          cout << a[i+k*3];
          if (i<2) cout << " , ";
        }
        cout << "> " << endl;
      }
      cout << "\tf = " << s->f << endl;
      for (int i=0;i<n;i++) {
        cout << "\tgradient[" << i << "] = "
             << gsl_vector_get(s->gradient,i) << endl;
      }
//  }
    gsl_multimin_fdfminimizer_free(s);
  }
*/
/*
  gsl_multimin_function my_func2;
  my_func2.n=n;
  my_func2.f=&lj;
  my_func2.params=reinterpret_cast<void*>(&atoms);
  for (int it2=0;it2<20;it2++) {
    for (int i=0;i<n;i++) {
      gsl_vector_set(xx,i,gsl_vector_get(x,i)*(1.+drand48()*1.e1));
//    gsl_vector_set(xx,i,gsl_vector_get(x,i));
    }
    gsl_vector *ss=gsl_vector_alloc(n); // step sizes
    gsl_vector_set_all(ss,1.e0);
    const gsl_multimin_fminimizer_type *T=
      gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s=gsl_multimin_fminimizer_alloc(T,n);
    gsl_multimin_fminimizer_set(s,&my_func2,xx,ss);
    int it=0;
//  see gsl_errno.h for status codes
    int status=GSL_CONTINUE;
    do {
//    Tracer t("iteration loop");
      it++;
      status=gsl_multimin_fminimizer_iterate(s);
//    cout << "\tminimizer status = " << status << endl;
//    if (status!=GSL_CONTINUE) break;
      double size=gsl_multimin_fminimizer_size(s);
//    status=gsl_multimin_test_size(size,1.e-2);
      status=GSL_CONTINUE;

//    for (int i=0;i<n;i++) {
//      cout << "\tx[" << i << "] = " << gsl_vector_get(s->x,i) << endl;
//    }
//    cout << "\tf = " << s->f << endl;
    } while (status==GSL_CONTINUE && it<itmax);
    cout << "\tstatus = " << status << endl;
//  if (status==GSL_SUCCESS) {
      cout << "\ttrial = " << it2 << endl;
      aFromX(atoms,s->x,a);
      for (int k=0;k<atoms;k++) {
        cout << "atom[" << k << "] = <";
        for (int i=0;i<3;i++) {
          cout << a[i+k*3];
          if (i<2) cout << " , ";
        }
        cout << "> " << endl;
      }
      cout << "\tf = " << s->fval << endl;
      gsl_vector *gggg=gsl_vector_alloc(n);
      ljgrad(s->x,&atoms,gggg);
      for (int i=0;i<n;i++) {
        cout << "\tgradient[" << i << "] = "
             << gsl_vector_get(gggg,i) << endl;
      }
      gsl_vector_free(gggg);
//  }
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
  }
  gsl_vector_free(x);
  gsl_vector_free(xx);
*/
  delete [] a; a=0;

  x=gsl_vector_alloc(1);
  gsl_vector_set(x,0,pow(2.,1./6.));
  a=new double[6];
  aFromX(2,x,a);
  drand48();
  for (int atoms=3;atoms<=5;atoms++) {
    Tracer t("atoms loop");
//  for (int count=0;count<20;count++) {
//    cout << "\n\tcount = " << count << endl;
      double *aa=new double[3*atoms];
      int atoms1=atoms-1;
      for (int i=0;i<3*atoms1;i++) aa[i]=a[i];
      for (int i=atoms*3-3;i<atoms*3;i++) aa[i]=a[i-3];
      if (atoms==3) aa[7]+=10.;
      else if (atoms==4) aa[atoms*3-1]+=10.;
      else aa[atoms*3-1]=-aa[atoms*3-4];
//    for (int k=0;k<atoms;k++) {
//      cout << "\tatom[" << k << "] = <";
//      for (int i=0;i<3;i++) {
//        cout << aa[i+k*3];
//        if (i<2) cout << " , ";
//      }
//      cout << "> " << endl;
//    }
      int n=max(1,3*atoms-6);
      gsl_vector *xx=gsl_vector_alloc(n);
      gsl_vector_set(xx,0,aa[3]);
      gsl_vector_set(xx,1,aa[6]);
      gsl_vector_set(xx,2,aa[7]);
      for (int i=3;i<n;i++) gsl_vector_set(xx,i,aa[i+6]);
//    cout << "\tenergy = " << lj(xx,&atoms) << endl;

/*
      gsl_multimin_function_fdf my_func;
      my_func.n=n;
      my_func.f=&lj;
      my_func.df=&ljgrad;
      my_func.fdf=&ljfdf;
      my_func.params=reinterpret_cast<void*>(&atoms);
      const gsl_multimin_fdfminimizer_type *T=
        gsl_multimin_fdfminimizer_vector_bfgs2;
      gsl_multimin_fdfminimizer *s=gsl_multimin_fdfminimizer_alloc(T,n);
      gsl_multimin_fdfminimizer_set(s,&my_func,xx,step_size,tol);
*/
//
      gsl_vector *ss=gsl_vector_alloc(n); // step sizes
      gsl_vector_set_all(ss,1.e0);
      const gsl_multimin_fminimizer_type *T=
        gsl_multimin_fminimizer_nmsimplex2;
      gsl_multimin_fminimizer *s=gsl_multimin_fminimizer_alloc(T,n);
      gsl_multimin_function my_func3;
      my_func3.n=n;
      my_func3.f=&lj;
      my_func3.params=reinterpret_cast<void*>(&atoms);
      gsl_multimin_fminimizer_set(s,&my_func3,xx,ss);
//

      int it=0;
      int status=GSL_CONTINUE;
      do {
        it++;
//      status=gsl_multimin_fdfminimizer_iterate(s);
        status=gsl_multimin_fminimizer_iterate(s);
      } while (/*status==GSL_CONTINUE && */ it<itmax);
      cout << "\tstatus = " << status << endl;
//    if (status==GSL_SUCCESS) {
        aFromX(atoms,s->x,aa);
        for (int k=0;k<atoms;k++) {
          cout << "atom[" << k << "] = <";
          for (int i=0;i<3;i++) {
            cout << aa[i+k*3];
            if (i<2) cout << " , ";
          }
          cout << "> " << endl;
        }
//      cout << "\tf = " << s->f << endl;
        cout << "\tf = " << s->fval << endl;
        gsl_vector *gggg=gsl_vector_alloc(n);
        ljgrad(s->x,&atoms,gggg);
        for (int i=0;i<n;i++) {
          cout << "\tgradient[" << i << "] = "
               << gsl_vector_get(gggg,i) << endl;
        }
        gsl_vector_free(gggg);
//    }
/*
      gsl_multimin_fminimizer_free(s);
*/
//
      gsl_vector_free(ss);
      gsl_multimin_fminimizer_free(s);
//
      gsl_vector_free(x);
      x=xx;
//    if (count < 19) {
//      delete [] aa; aa=0;
//    } else
      a=aa;
//  }
  }
  delete [] a; a=0;
  gsl_vector_free(x);

/*
  for (int i=0;i<n;i++) gsl_vector_set(x,i,gsl_vector_get(s->x,i));
  double energy=lj(x,&atoms);
  gsl_vector *g=gsl_vector_alloc(n);
  gsl_vector *gg=gsl_vector_alloc(n);
  ljgrad(x,&atoms,g);
  double ljfdf_energy=0.;
  ljfdf(x,&atoms,&ljfdf_energy,gg);
  cout << "\tenergy = " << energy << " " << ljfdf_energy << endl;
  double dx=sqrt(numeric_limits<double>::epsilon());
  for (int i=0;i<n;i++) {
    double xi=gsl_vector_get(x,i);
    gsl_vector_set(x,i,xi+dx);
    double energyn=lj(x,&atoms);
    double deriv=(energyn-energy)/dx;
    cout << "\tg[" << i << "] = " << gsl_vector_get(g,i) << " "
         << gsl_vector_get(gg,i) << " " << deriv << " "
         << ( gsl_vector_get(g,i) - deriv ) / gsl_vector_get(g,i) << endl;
    gsl_vector_set(x,i,xi);
  }
  gsl_vector_free(g);
  gsl_vector_free(gg);
*/
}
