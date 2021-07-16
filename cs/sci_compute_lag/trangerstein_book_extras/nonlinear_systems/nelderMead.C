#include "Vector.H"

void selectionSort(int n,double *fx,Vector<double> **x) {
  for (int j = 0; j < n; j++) {
    int k = j;
    for (int i = j+1; i <= n; i++) {
      if (fx[i] < fx[k]) k = i;
    }
    if ( k != j ) {
      double t=fx[j];
      fx[j]=fx[k];
      fx[k]=t;
      Vector<double> *xt=x[j];
      x[j]=x[k];
      x[k]=xt;
   }
 }
}

void nelderMead(int maxit,double tol,double (*f)(const Vector<double>&),
Vector<double> **x) {
  int count=0;
  int n=x[0]->size();
  double *fx=OPERATOR_NEW_BRACKET(double,n+1);
  for (int i=0;i<=n;i++) {
    fx[i]=f(*x[i]); count++;
  }
  selectionSort(n,fx,x);
  while (count < maxit && fx[n]-fx[0]>tol) {
    Vector<double> xave(n,0.);
    for (int i=0;i<n;i++) xave+= *x[i];
    xave /= static_cast<double>(n);
    Vector<double> x_reflect= ( xave * 2. ) - *x[n];
    double f_reflect=f(x_reflect); count++;

    if (f_reflect < fx[0]) {
//    expand
      Vector<double> x_expand=(xave*3)-(*x[n]*2);
      double f_expand=f(x_expand); count++
      if (f_expand < f_reflect) {
        *x[n]=x_expand;
        fx[n]=f_expand;
      } else {
        *x[n]=x_reflect;
        fx[n]=f_reflect;
      }
    } else if (f_reflect < fx[n-1]) {
//    reflect
      *x[n]=x_reflect;
        fx[n]=f_reflect;
    } else {
      bool shrink=false;
      if (f_reflect < fx[n]) {
//      outside contraction
        Vector<double> x_oc=(xave*1.5)-(*x[n]*.5);
        double f_oc=f(x_oc); count++
        if (f_oc<=f_reflect) {
          *x[n]=x_oc;
          fx[n]=f_oc;
        } else shrink=true;
      } else {
//      inside contraction
        Vector<double x_ic=(xave+*x[n])*.5;
        double f_ic=f(x_ic); count++
        if (f_ic<fx[n]) {
          *x[n]=x_ic;
          fx[n]=f_ic;
        } else shrink=true;
      }
      if (shrink) {
        for (int i=1;i<=n;i++) {
          *x[i]=(*x[0]*1.5)-(*x[i]*.5);
          fx[i]=f(*x[i]); count++;
        }
      }
    }
    selectionSort(n,fx,x);
  }
}
