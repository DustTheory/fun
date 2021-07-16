// "$Header:$"
//----------------------------  polynomial.cc  -----------------------
//      $Id: polynomial.cc,v 1.35 2003/05/12 21:37:34 wolf Exp $   
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------   polynomial.cc  ----------------------
//
//modified from deal.II/base/source/polynomial.cc
//  by John Trangenstein, August 2009
//**********************************************************************
// Copyright 2009 John A. Trangenstein
//
// This software is made available for research and instructional use 
// only. 
// You may copy and use this software without charge for these 
// non-commercial purposes, provided that the copyright notice and 
// associated text is reproduced on all copies.  
// For all other uses (including distribution of modified versions), 
// please contact the author at
//   John A. Trangenstein
//   Department of Mathematics
//   Duke University
//   Durham, NC 27708-0320
//   USA
// or
//   johnt@math.duke.edu
// 
// This software is made available "as is" without any assurance that it
// is completely correct, or that it will work for your purposes.  
// Use the software at your own risk.
//**********************************************************************

#include <algorithm>
#include <cmath>
#include <Const.H>
#include <MemoryDebugger.H>
#include <Polynomial.H>
#include <Tracer.H>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Polynomial::printOn(ostream &os) const {
  os << "Polynomial:" << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Monomial::values(REAL x,NumPtr<REAL> &v) const {
  CHECK_SAME(order+1,v.getNumber());
  v[0]=1.;
  for (int i=1;i<=order;i++) v[i]=x*v[i-1];
}

void Monomial::slopes(REAL x,NumPtr<REAL> &s) const {
  CHECK_SAME(order+1,s.getNumber());
  s[0]=0.;
  REAL xi=1.;
  for (int i=1;i<=order;i++) {
    s[i]=static_cast<REAL>(i)*xi;
    xi*=x;
  }
}

void Monomial::slope2s(REAL x,NumPtr<REAL> &s2) const {
  CHECK_SAME(order+1,s2.getNumber());
  s2[0]=0.;
  if (order>=1) s2[1]=0.;
  REAL xi=1.;
  for (int i=2;i<=order;i++) {
    s2[i]=static_cast<REAL>(i*(i-1))*xi;
    xi*=x;
  }
}

//values[j] = j'th derivative
void Monomial::derivatives(REAL x,int j,NumPtr<REAL> &values) const {
  int number_derivs=values.getNumber();
  if (number_derivs<=0) return;
  CHECK_TEST(number_derivs<=3);
  REAL xp=1.;
  for (int i=0;i<j-2;i++) xp*=x;
  if (number_derivs>=3) values[2]=static_cast<REAL>(j*(j-1))*xp;
  xp*=x;
  if (number_derivs>=2) values[1]=static_cast<REAL>(j)*xp;
  values[0]=xp*x;
}

void Monomial::printOn(ostream &os) const {
  os << "Monomial: order = " << order << endl;
  Polynomial::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C0LagrangePolynomial::C0LagrangePolynomial(int ord) : order(ord),
den(ord+1),zeros(ord+1) {
  if (order>0) {
    REAL ro=1./static_cast<REAL>(order);
    for (int i=0;i<=order;i++) zeros[i]=static_cast<REAL>(i)*ro;
  } else {
    zeros[0]=0.;
  }
  for (int j=0;j<=order;j++) {
    REAL denj=1.;
    for (int i=0;i<=order;i++) {
      if (i!=j) denj*=zeros[j]-zeros[i];
    }
    den[j]=denj;
  }
}

void C0LagrangePolynomial::values(REAL x,NumPtr<REAL> &v) const {
  for (int j=0;j<=order;j++) {
    REAL num=1.;
    for (int i=0;i<=order;i++) {
      if (i!=j) num*=x-zeros[i];
    }
    v[j]=num/den[j];
  }
}

void C0LagrangePolynomial::slopes(REAL x,NumPtr<REAL> &s) const {
  for (int j=0;j<=order;j++) {
    REAL sum=0.;
    for (int k=0;k<=order;k++) {
      if (k==j) continue;
      REAL prod=1.;
      for (int i=0;i<=order;i++) {
        if (i!=j && i!=k) prod*=x-zeros[i];
      }
      sum+=prod;
    }
    s[j]=sum/den[j];
  }
}

void C0LagrangePolynomial::slope2s(REAL x,NumPtr<REAL> &s2) const {
  for (int j=0;j<=order;j++) {
    REAL sum=0.;
    for (int k=0;k<=order;k++) {
      if (k==j) continue;
      for (int m=0;m<=order;m++) {
        if (m==j || m==k) continue;
        REAL prod=1.;
        for (int i=0;i<=order;i++) {
          if (i!=j && i!=k && i!=m) prod*=x-zeros[i];
        }
        sum+=prod;
      }
    }
    s2[j]=sum/den[j];
  }
}

//values[j] = j'th derivative
void C0LagrangePolynomial::derivatives(REAL x,int j,
NumPtr<REAL> &values) const {
  int number_derivs=values.getNumber();
  if (number_derivs<=0) return;
  CHECK_TEST(number_derivs<=3);
  REAL num=1.;
  REAL sum1=0.;
  REAL sum2=0.;
  for (int k=0;k<=order;k++) {
    if (k==j) continue;
    num*=x-zeros[k];
    if (number_derivs<2) continue;
    REAL prod1=1.;
    for (int m=0;m<=order;m++) {
      if (m==j || m==k) continue;
      prod1*=x-zeros[m];
      if (number_derivs<3) continue;
      REAL prod2=1.;
      for (int i=0;i<=order;i++) {
        if (i!=j && i!=k && i!=m) prod2*=x-zeros[i];
      }
      sum2+=prod2;
    }
    sum1+=prod1;
  }
  values[0]=num/den[j];
  if (number_derivs>=2) values[1]=sum1/den[j];
  if (number_derivs>=3) values[2]=sum2/den[j];
}

void C0LagrangePolynomial::printOn(ostream &os) const {
  os << "C0LagrangePolynomial: order = " << order << endl;
  Polynomial::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C1LagrangePolynomial::C1LagrangePolynomial(int ord) {
  n=max(1,(ord-1)/2);
  int np1=n+1;
  order=2*n+1;
  zeros.allocate(np1);
  den.allocate(np1);
  a.allocate(np1);
  b.allocate(np1);
  REAL ro=1./static_cast<REAL>(n);
  for (int i=0;i<=n;i++) zeros[i]=static_cast<REAL>(i)*ro;
  for (int j=0;j<=n;j++) {
    REAL sum=0.;
    for (int i=0;i<=n;i++) {
      if (i!=j) sum+=1./(zeros[j]-zeros[i]);
    }
    a[j]=-2.*sum;
    b[j]=1.-a[j]*zeros[j];
    REAL denj=1.;
    for (int i=0;i<=n;i++) {
      if (i!=j) denj*=pow(zeros[j]-zeros[i],2);
    }
    den[j]=denj;
  }
}

void C1LagrangePolynomial::values(REAL x,NumPtr<REAL> &v) const {
  CHECK_SAME(order+1,v.getNumber());
  int np1=n+1;
  for (int j=0;j<=n;j++) {
    REAL num=1.;
    for (int i=0;i<=n;i++) {
      if (i!=j) num*=pow(x-zeros[i],2);
    }
    REAL quot=num/den[j];
    v[j]=(a[j]*x+b[j])*quot;
    v[j+np1]=(x-zeros[j])*quot;
  }
}

void C1LagrangePolynomial::slopes(REAL x,NumPtr<REAL> &s) const {
  CHECK_SAME(order+1,s.getNumber());
  int np1=n+1;
  for (int j=0;j<=n;j++) {
    REAL num=1.;
    REAL sum=0.;
    for (int k=0;k<=n;k++) {
      if (k==j) continue;
      num*=pow(x-zeros[k],2);
      REAL prod=2.*(x-zeros[k]);
      for (int i=0;i<=n;i++) {
        if (i!=j && i!=k) prod*=pow(x-zeros[i],2);
      }
      sum+=prod;
    }
    s[j]=(a[j]*num+(a[j]*x+b[j])*sum)/den[j];
    s[j+np1]=(num+(x-zeros[j])*sum)/den[j];
  }
}

void C1LagrangePolynomial::slope2s(REAL x,NumPtr<REAL> &s2) const {
  CHECK_SAME(order+1,s2.getNumber());
  int np1=n+1;
  for (int j=0;j<=n;j++) {
    REAL sum1=0.;
    REAL sum2=0.;
    REAL sum3=0.;
    for (int k=0;k<=n;k++) {
      if (k==j) continue;
      REAL num=1.;
      REAL sum=0.;
      for (int m=0;m<=n;m++) {
        if (m==j || m==k) continue;
        num*=pow(x-zeros[m],2);
        REAL prod=1.;
        for (int i=0;i<=n;i++) {
          if (i!=j && i!=k && i!=m) prod*=pow(x-zeros[i],2);
        }
        sum+=(x-zeros[m])*prod;
      }
      sum1+=(x-zeros[k])*num;
      sum2+=num;
      sum3+=(x-zeros[k])*sum;
    }
    sum1*=4.;
    sum2*=2.;
    sum3*=4.;
    s2[j]=(a[j]*sum1+(a[j]*x+b[j])*(sum2+sum3))/den[j];
    s2[j+np1]=(sum1+(x-zeros[j])*(sum2+sum3))/den[j];
  }
}

//values[j] = j'th derivative
void C1LagrangePolynomial::derivatives(REAL x,int J,
NumPtr<REAL> &values) const {
  int number_derivs=values.getNumber();
  if (number_derivs<=0) return;
  CHECK_TEST(number_derivs<=3);

  int j=J%(n+1);
  REAL num=1.;
  REAL sum1=0.;
  REAL sum2=0.;
  REAL sum3=0.;
  for (int k=0;k<=n;k++) {
    if (k==j) continue;
    num*=pow(x-zeros[k],2);
    if (number_derivs<2) continue;
    REAL prod1=1.;
    REAL sum=0.;
    for (int m=0;m<=n;m++) {
      if (m==j || m==k) continue;
      prod1*=pow(x-zeros[m],2);
      if (number_derivs<3) continue;
      REAL prod2=1.;
      for (int i=0;i<=n;i++) {
        if (i!=j && i!=k && i!=m) prod2*=pow(x-zeros[i],2);
      }
      sum+=(x-zeros[m])*prod2;
    }
    sum1+=(x-zeros[k])*prod1;
    sum2+=prod1;
    sum3+=(x-zeros[k])*sum;
  }
  sum1*=2.;
  sum2*=2.;
  sum3*=4.;
  if (J<=n) {
    REAL sum=0.;
    for (int i=0;i<=n;i++) {
      if (i!=j) sum+=1./(zeros[j]-zeros[i]);
    }
    REAL axb=a[j]*x+b[j];
    values[0]=axb*num/den[j];
    if (number_derivs>=2) values[1]=(a[j]*num+axb*sum1)/den[j];
    if (number_derivs>=3) {
      values[2]=(2.*a[j]*sum1+axb*(sum2+sum3))/den[j];
    }
  } else {
    REAL xz=x-zeros[j];
    values[0]=xz*num/den[j];
    if (number_derivs>=2) values[1]=(num+xz*sum1)/den[j];
    if (number_derivs>=3) values[2]=(2.*sum1+xz*(sum2+sum3))/den[j];
  }
}

void C1LagrangePolynomial::printOn(ostream &os) const {
  os << "C1LagrangePolynomial: order = " << order << endl;
  Polynomial::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C2LagrangePolynomial::C2LagrangePolynomial(int ord) {
  n=max(1,(ord-2)/3);
  int np1=n+1;
  order=3*n+2;
  zeros.allocate(np1);
  a0.allocate(np1);
  b0.allocate(np1);
  c0.allocate(np1);
  a1.allocate(np1);
  b1.allocate(np1);
  den.allocate(np1);
  REAL ro=1./static_cast<REAL>(n);
  for (int i=0;i<=n;i++) zeros[i]=static_cast<REAL>(i)*ro;
  for (int j=0;j<=n;j++) {
    REAL sum1=0.;
    REAL sum2=0.;
    REAL sum3=0.;
    for (int i=0;i<=n;i++) {
      if (i==j) continue;
      REAL sum=0.;
      for (int m=0;m<=n;m++) {
        if (m!=i && m!=j) sum+=1./(zeros[j]-zeros[m]);
      }
      REAL dz=1./(zeros[j]-zeros[i]);
      sum1+=dz;
      sum2+=pow(dz,2);
      sum3+=sum*dz;
    }
    a0[j]=4.5*(pow(sum1,2)-sum3)-3.*sum2;
    b0[j]=-3.*sum1-2.*a0[j]*zeros[j];
    c0[j]=1.-zeros[j]*(b0[j]+a0[j]*zeros[j]);
    a1[j]=-3.*sum1;
    b1[j]=1.-a1[j]*zeros[j];
    REAL denj=1.;
    for (int i=0;i<=n;i++) {
      if (i!=j) denj*=pow(zeros[j]-zeros[i],3);
    }
    den[j]=denj;
  }
}

void C2LagrangePolynomial::values(REAL x,NumPtr<REAL> &v) const {
  CHECK_SAME(order+1,v.getNumber());
  int np1=n+1;
  int np12=2*np1;
  for (int j=0;j<=n;j++) {
    REAL num=1.;
    for (int i=0;i<=n;i++) {
      if (i!=j) num*=pow(x-zeros[i],3);
    }
    REAL quot=num/den[j];
    v[j]=((a0[j]*x+b0[j])*x+c0[j])*quot;
    REAL xz=x-zeros[j];
    v[j+np1]=(a1[j]*x+b1[j])*xz*quot;
    v[j+np12]=0.5*xz*xz*quot;
  }
}

void C2LagrangePolynomial::slopes(REAL x,NumPtr<REAL> &s) const {
  CHECK_SAME(order+1,s.getNumber());
  int np1=n+1;
  int np12=2*np1;
  for (int j=0;j<=n;j++) {
    REAL num=1.;
    REAL sum=0.;
    for (int k=0;k<=n;k++) {
      if (k==j) continue;
      REAL xz=x-zeros[k];
      REAL prod=xz*xz;
      num*=prod*xz;
      for (int i=0;i<=n;i++) {
        if (i!=j && i!=k) prod*=pow(x-zeros[i],3);
      }
      sum+=prod;
    }
    num/=den[j];
    sum*=3./den[j];
    s[j]=(2.*a0[j]*x+b0[j])*num+((a0[j]*x+b0[j])*x+c0[j])*sum;
    REAL axb1=a1[j]*x+b1[j];
    REAL xz=x-zeros[j];
    s[j+np1]=(a1[j]*xz+axb1)*num+axb1*xz*sum;
    s[j+np12]=xz*(num+0.5*xz*sum);
  }
}

void C2LagrangePolynomial::slope2s(REAL x,NumPtr<REAL> &s2) const {
  CHECK_SAME(order+1,s2.getNumber());
  int np1=n+1;
  int np12=2*np1;
  for (int j=0;j<=n;j++) {
    REAL num=1.;
    REAL sum1=0.;
    REAL sum2=0.;
    REAL sum3=0.;
    for (int k=0;k<=n;k++) {
      if (k==j) continue;
      REAL num12=1.;
      REAL sum=0.;
      for (int m=0;m<=n;m++) {
        if (m==j || m==k) continue;
        REAL prod=1.;
        for (int i=0;i<=n;i++) {
          if (i!=j && i!=k && i!=m) prod*=pow(x-zeros[i],3);
        }
        REAL xzm=x-zeros[m];
        REAL xzm2=xzm*xzm;
        num12*=xzm2*xzm;
        sum+=xzm2*prod;
      }
      REAL xz=x-zeros[k];
      REAL xz2=xz*xz;
      num*=xz2*xz;
      sum1+=xz2*num12;
      sum2+=xz*num12;
      sum3+=xz2*sum;
    }
    sum1*=3.;
    sum2*=6.;
    sum3*=9.;
    s2[j]=(2.*(a0[j]*num+(2.*a0[j]*x+b0[j])*sum1)
          +((a0[j]*x+b0[j])*x+c0[j])*(sum2+sum3))/den[j];
    REAL xz=x-zeros[j];
    REAL axb1=a1[j]*x+b1[j];
    s2[j+np1]=(2.*a1[j]*num+2.*(a1[j]*xz+axb1)*sum1
              + axb1*xz*(sum2+sum3))/den[j];
    s2[j+np12]=(num+xz*(2.*sum1+0.5*xz*(sum2+sum3)))/den[j];
  }
}

//values[j] = j'th derivative
void C2LagrangePolynomial::derivatives(REAL x,int J,
NumPtr<REAL> &values) const {
  int number_derivs=values.getNumber();
  if (number_derivs<=0) return;
  CHECK_TEST(number_derivs<=3);

  int np1=n+1;
  int j=J%(np1);
  REAL num=1.;
  REAL sum1=0.;
  REAL sum2=0.;
  REAL sum3=0.;
  for (int k=0;k<=n;k++) {
    if (k==j) continue;
    num*=pow(x-zeros[k],3);
    if (number_derivs<2) continue;
    REAL prod1=1.;
    REAL sum=0.;
    for (int m=0;m<=n;m++) {
      if (m==j || m==k) continue;
      REAL xzm=x-zeros[m];
      prod1*=pow(xzm,3);
      if (number_derivs<3) continue;
      REAL prod2=1.;
      for (int i=0;i<=n;i++) {
        if (i!=j && i!=k && i!=m) prod2*=pow(x-zeros[i],3);
      }
      sum+=xzm*xzm*prod2;
    }
    REAL xzk=x-zeros[k];
    REAL xzk2=xzk*xzk;
    sum1+=xzk2*prod1;
    sum2+=xzk*prod1;
    sum3+=xzk2*sum;
  }
  sum1*=3.;
  sum2*=6.;
  sum3*=9.;
  if (J<np1) {
    REAL ax=a0[j]*x;
    REAL abcx=ax+b0[j];
    REAL abx=abcx+ax;
    abcx=c0[j]+x*abcx;
    values[0]=abcx*num/den[j];
    if (number_derivs>=2) values[1]=(abx*num+abcx*sum1)/den[j];
    if (number_derivs>=3) {
      values[2]=(2.*a0[j]*num+2.*abx*sum1+abcx*(sum2+sum3))/den[j];
    }
  } else if (j<np1+np1) {
    REAL xz=x-zeros[j];
    REAL axb=a1[j]*x+b1[j];
    values[0]=axb*xz*num/den[j];
    if (number_derivs>=2) {
      values[1]=((a1[j]*xz+axb)*num+axb*xz*sum1)/den[j];
    }
    if (number_derivs>=3) {
      values[2]=(2.*(a1[j]*num+(a1[j]*xz+axb)*sum1)
                +axb*xz*(sum2+sum3))/den[j];
    }
  } else {
    REAL xz=x-zeros[j];
    values[0]=0.5*num/den[j];
    if (number_derivs>=2) values[1]=xz*(num+0.5*xz*sum1)/den[j];
    if (number_derivs>=3) {
      values[2]=(num+xz*(2.*sum1+0.5*xz*(sum2+sum3)))/den[j];
    }
  }
}

void C2LagrangePolynomial::printOn(ostream &os) const {
  os << "C2LagrangePolynomial: order = " << order << endl;
  Polynomial::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//on input, assumes x in [0,1]
void LegendrePolynomial::values(REAL x,NumPtr<REAL> &v) const {
  int order=v.getNumber()-1;
  CHECK_TEST(order>=0);
  v[0]=1.;
  if (order<1) return;

  x=2.*x-1.;
  REAL pm1=x;
  v[1]=pm1;
  if (order<2) return;

  REAL pm2=1.;
  REAL r2im1=1.;
  REAL ri=1.;
  for (int i=2;i<=order;i++) {
    r2im1+=2.;
    ri+=1.;
    v[i]=(r2im1*x*pm1-(ri-1.)*pm2)/ri;
    pm2=pm1;
    pm1=v[i];
  }
}

void LegendrePolynomial::slopes(REAL x,NumPtr<REAL> &s) const {
  int order=s.getNumber()-1;
  CHECK_TEST(order>=0);
  s[0]=0.;
  if (order<1) return;

  REAL dpm1=1.;
  s[1]=2.*dpm1;
  if (order<2) return;

  REAL dpm2=0.;
  x=2.*x-1.;
  REAL pm2=1.;
  REAL pm1=x;
  REAL r2im1=1.;
  REAL ri=1.;
  for (int i=2;i<=order;i++) {
    r2im1+=2.;
    ri+=1.;
    REAL p=(r2im1*x*pm1-(ri-1.)*pm2)/ri;
    REAL dp=r2im1*pm1+dpm2;
    s[i]=2.*dp;
    pm2=pm1;
    pm1=p;
    dpm2=dpm1;
    dpm1=dp;
  }
}

void LegendrePolynomial::slope2s(REAL x,NumPtr<REAL> &s2) const {
  int order=s2.getNumber()-1;
  CHECK_TEST(order>=0);
  s2[0]=0.;
  if (order<1) return;

  s2[1]=0.;
  if (order<2) return;

  x=2.*x-1.;
  REAL pm2=1.;
  REAL pm1=x;
  REAL r2im1=1.;
  REAL ri=1.;
  REAL dpm2=0.;
  REAL dpm1=1.;
  REAL ddpm1=0.;
  REAL ddpm2=0.;
  for (int i=2;i<=order;i++) {
    r2im1+=2.;
    ri+=1.;
    REAL p=(r2im1*x*pm1-(ri-1.)*pm2)/ri;
    REAL dp=(r2im1*(pm1+x*dpm1)-(ri-1.)*dpm2)/ri;
    REAL ddp=(r2im1*(2.*dpm1+x*ddpm1)-(ri-1.)*ddpm2)/ri;
    s2[i]=4.*ddp;
    pm2=pm1;
    pm1=p;
    dpm2=dpm1;
    dpm1=dp;
    ddpm2=ddpm1;
    ddpm1=ddp;
  }
}

void LegendrePolynomial::derivatives(REAL x,int order,
NumPtr<REAL> &values) const {
  int number_derivs=values.getNumber();
  if (number_derivs<=0) return;
  CHECK_TEST(number_derivs<=3);
  if (order==0) {
    values[0]=1.;
    if (number_derivs>1) values[1]=0.;
    if (number_derivs>2) values[2]=0.;
    return;
  }
  x=2.*x-1.;
  values[0]=x;
  if (number_derivs>1) values[1]=1.;
  if (number_derivs>2) values[2]=0.;
  REAL pm2=1.;
  REAL dpm2=0.;
  REAL ddpm2=0.;
  REAL r2im1=1.;
  REAL ri=1.;
  for (int i=2;i<=order;i++) {
    r2im1+=2.;
    ri+=1.;
    REAL p=(r2im1*x*values[0]-(ri-1.)*pm2)/ri;
    if (number_derivs>1) {
      REAL dp=(r2im1*(values[0]+x*values[1])-(ri-1.)*dpm2)/ri;
      if (number_derivs>2) {
        REAL ddp=(r2im1*(2.*values[1]+x*values[2])-(ri-1.)*ddpm2)/ri;
        ddpm2=values[2];
        values[2]=ddp;
      }
      dpm2=values[1];
      values[1]=dp;
    }
    pm2=values[0];
    values[0]=p;
  }
  if (number_derivs>1) values[1]*=2.;
  if (number_derivs>2) values[2]*=4.;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void LegendrePolynomial::printOn(ostream &os) const {
  os << "LegendrePolynomial:" << endl;
  Polynomial::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HierarchicalPolynomial::values(REAL x,NumPtr<REAL> &v) const {
  int order=v.getNumber()-1;
  CHECK_TEST(order>=1);
  x=2.*x-1.;
  v[0]=0.5*(1.-x);
  v[1]=0.5*(1.+x);
  if (order<2) return;
  REAL pm2=1.;
  REAL pm1=x;
  REAL r2im1=1.;
  REAL ri=1.;
  for (int i=2;i<=order;i++) {
    r2im1+=2.;
    ri+=1.;
    REAL p=(r2im1*x*pm1-(ri-1.)*pm2)/ri;
    v[i]=0.5*(p-pm2)/sqrt(r2im1);
    pm2=pm1;
    pm1=p;
  }
}

void HierarchicalPolynomial::slopes(REAL x,NumPtr<REAL> &s) const {
  int order=s.getNumber()-1;
  CHECK_TEST(order>=1);
  s[0]=-1.;
  s[1]=1.;
  if (order<2) return;
  x=2.*x-1.;
  REAL pm2=1.;
  REAL pm1=x;
  REAL r2im1=1.;
  REAL ri=1.;
  for (int i=2;i<=order;i++) {
    r2im1+=2.;
    ri+=1.;
    s[i]=sqrt(r2im1)*pm1;
    REAL p=(r2im1*x*pm1-(ri-1.)*pm2)/ri;
    pm2=pm1;
    pm1=p;
  }
}

void HierarchicalPolynomial::slope2s(REAL x,NumPtr<REAL> &s) const {
  int order=s.getNumber()-1;
  CHECK_TEST(order>=1);
  s[0]=0.;
  s[1]=0.;
  if (order<2) return;
  x=2.*x-1.;
  REAL pm2=1.;
  REAL pm1=x;
  REAL r2im1=1.;
  REAL ri=1.;
  REAL dpm2=0.;
  REAL dpm1=1.;
  for (int i=2;i<=order;i++) {
    r2im1+=2.;
    ri+=1.;
    REAL p=(r2im1*x*pm1-(ri-1.)*pm2)/ri;
    REAL dp=r2im1*pm1+dpm2;
    s[i]=2.*dpm1*sqrt(r2im1);
    pm2=pm1;
    pm1=p;
    dpm2=dpm1;
    dpm1=dp;
  }
}

void HierarchicalPolynomial::derivatives(REAL x,int order,
NumPtr<REAL> &values) const {
  int number_derivs=values.getNumber();
  if (number_derivs<=0) return;
  CHECK_TEST(number_derivs<=3);
  x=2.*x-1.;
  if (order==0) {
    values[0]=0.5*(1.-x);
    if (number_derivs>1) values[1]=-1.;
    if (number_derivs>2) values[2]=0.;
    return;
  }
  if (order==1) {
    values[0]=0.5*(1.+x);
    if (number_derivs>1) values[1]=1.;
    if (number_derivs>2) values[2]=0.;
    return;
  }
  REAL pm1=x;
  REAL dpm1=1.;
  REAL ddpm1=0.;
  REAL pm2=1.;
  REAL dpm2=0.;
  REAL ddpm2=0.;
  REAL r2im1=1.;
  REAL ri=1.;
  for (int i=2;i<=order;i++) {
    r2im1+=2.;
    ri+=1.;
    REAL p=(r2im1*x*pm1-(ri-1.)*pm2)/ri;
    REAL dp=UNDEFINED;
    REAL ddp=UNDEFINED;
    if (number_derivs>1) {
      dp=(r2im1*(pm1+x*dpm1)-(ri-1.)*dpm2)/ri;
      if (number_derivs>2) {
        ddp=(r2im1*(2.*dpm1+x*ddpm1)-(ri-1.)*ddpm2)/ri;
      }
    }
    if (i==order) {
      REAL root=sqrt(r2im1);
      values[0]=0.5*(p-pm2)/root;
      if (number_derivs>1) values[1]=(dp-dpm2)/root;
      if (number_derivs>2) values[2]=(ddp-ddpm2)*2./root;
      return;
    }
    if (number_derivs>2) {
      ddpm2=ddpm1;
      ddpm1=ddp;
    }
    if (number_derivs>1) {
      dpm2=dpm1;
      dpm1=dp;
    }
    pm2=pm1;
    pm1=p;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HierarchicalPolynomial::printOn(ostream &os) const {
  os << "HierarchicalPolynomial:" << endl;
  Polynomial::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HierarchicalSidePolynomial::values(REAL x,NumPtr<REAL> &v) const {
  int order=v.getNumber()-1;
  CHECK_TEST(order>=0);
  v[0]=-sqrt(3.);
  if (order<1) return;
  x=2.*x-1.;
  v[1]=-sqrt(5.)*x;
  if (order<2) return;
  for (int i=2;i<=order;i++) {
    v[i]=(x*v[i-1]*sqrt(static_cast<REAL>(3+i*(8+4*i)))
      -v[i-2]*static_cast<REAL>(i-1)
          *sqrt(static_cast<REAL>(2*i+3)/static_cast<REAL>(2*i-1)))
      /static_cast<REAL>(i+2);
  }
}

void HierarchicalSidePolynomial::slopes(REAL x,NumPtr<REAL> &s) const {
  int order=s.getNumber()-1;
  CHECK_TEST(order>=0);
  s[0]=0.;
  if (order<1) return;
  s[1]=-sqrt(20.);
  if (order<2) return;
  x=2.*x-1.;
  REAL pm2=-sqrt(3.);
  REAL pm1=-sqrt(5.)*x;
  for (int i=2;i<=order;i++) {
    REAL a=sqrt(static_cast<REAL>(3+i*(8+4*i)));
    REAL b=static_cast<REAL>(i-1)
      *sqrt(static_cast<REAL>(2*i+3)/static_cast<REAL>(2*i-1));
    REAL c=1./static_cast<REAL>(i+2);
    REAL p=(x*pm1*a-pm2*b)*c;
    s[i]=(a*(2.*pm1+x*s[i-1])-b*s[i-2])*c;
    pm2=pm1;
    pm1=p;
  }
}

void HierarchicalSidePolynomial::slope2s(REAL x,NumPtr<REAL> &s) const {
  int order=s.getNumber()-1;
  CHECK_TEST(order>=0);
  s[0]=0.;
  if (order<1) return;
  s[1]=0.;
  if (order<2) return;
  x=2.*x-1.;
  REAL dpm2=0.;
  REAL dpm1=-sqrt(20.);
  REAL pm2=-sqrt(3.);
  REAL pm1=-sqrt(5.)*x;
  for (int i=2;i<=order;i++) {
    REAL a=sqrt(static_cast<REAL>(3+i*(8+4*i)));
    REAL b=static_cast<REAL>(i-1)
      *sqrt(static_cast<REAL>(2*i+3)/static_cast<REAL>(2*i-1));
    REAL c=1./static_cast<REAL>(i+2);
    REAL p=(x*pm1*a-pm2*b)*c;
    REAL dp=(a*(2.*pm1+x*dpm1)-b*dpm2)*c;
    s[i]=(a*(4.*dpm1+x*s[i-1])-b*s[i-2])*c;
    pm2=pm1;
    pm1=p;
    dpm2=dpm1;
    dpm1=dp;
  }
}

void HierarchicalSidePolynomial::derivatives(REAL x,int order,
NumPtr<REAL> &values) const {
  int number_derivs=values.getNumber();
  if (number_derivs<=0) return;
  CHECK_TEST(number_derivs<=3);
  x=2.*x-1.;
  if (order==0) {
    values[0]=-sqrt(3.);
    if (number_derivs>1) values[1]=0.;
    if (number_derivs>2) values[2]=0.;
    return;
  }
  if (order==1) {
    values[0]=-sqrt(5.)*x;
    if (number_derivs>1) values[1]=-sqrt(20.);;
    if (number_derivs>2) values[2]=0.;
    return;
  }
  REAL pm2=-sqrt(3.);
  REAL pm1=-sqrt(5.)*x;
  REAL dpm2=0.;
  REAL dpm1=-sqrt(20.);
  REAL ddpm2=0.;
  REAL ddpm1=0.;
  for (int i=2;i<=order;i++) {
    REAL a=sqrt(static_cast<REAL>(3+i*(8+4*i)));
    REAL b=static_cast<REAL>(i-1)
      *sqrt(static_cast<REAL>(2*i+3)/static_cast<REAL>(2*i-1));
    REAL c=1./static_cast<REAL>(i+2);
    REAL p=(x*pm1*a-pm2*b)*c;
    REAL dp=(a*(2.*pm1+x*dpm1)-b*dpm2)*c;
    REAL ddp=(a*(4.*dpm1+x*ddpm1)-b*ddpm2)*c;
    if (i==order) {
      values[0]=p;
      if (number_derivs>1) values[1]=dp;
      if (number_derivs>2) values[2]=ddp;
      return;
    }
    if (number_derivs>2) {
      ddpm2=ddpm1;
      ddpm1=ddp;
    }
    if (number_derivs>1) {
      dpm2=dpm1;
      dpm1=dp;
    }
    pm2=pm1;
    pm1=p;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HierarchicalSidePolynomial::printOn(ostream &os) const {
  os << "HierarchicalPolynomial:" << endl;
  Polynomial::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

template double __cmath_power<double>(double,unsigned);
#include <NumPtr.C>
INSTANTIATE_NUMPTR(REAL)
INSTANTIATE_NUMPTR(NumPtr<REAL>)
