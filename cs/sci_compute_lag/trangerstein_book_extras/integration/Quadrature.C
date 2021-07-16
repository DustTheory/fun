// "$Header:$"
//----------------------------  quadrature.cc  ---------------------------
//  $Id: quadrature.cc,v 1.45 2003/01/14 19:51:07 wolf Exp $
//  Version: $Name: Version-4-0-0 $
//
//  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//  This file is subject to QPL and may not be  distributed
//  without copyright and license information. Please refer
//  to the file deal.II/doc/license.html for the  text  and
//  further information on this license.
//
//----------------------------  quadrature.cc  ---------------------------
//------------------------  quadrature_lib.cc  ---------------------------
//    $Id: quadrature_lib.cc,v 1.32 2003/01/03 22:21:28 deal Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//------------------------  quadrature_lib.cc  ---------------------------
//
//modified from deal.II/base/source/quadrature.cc and quadrature_lib.cc
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

#include <Quadrature.H>
#include <Arch.H>
#include <NumPtr.H>
#include <Polynomial.H>
#include <Tensor.H>
#include <Tracer.H>
#include <Types.H>
#include <algorithm>
#include <cmath>
#include <iostream>
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <> Quadrature<0>::Quadrature(const Quadrature<-1> &,
const Quadrature<1> &) : number_points(0) {
  OBSOLETE("Quadrature impossible in zero dimensions");
}

template <> Quadrature<1>::Quadrature(const Quadrature<0> &,
const Quadrature<1> &) : number_points(0) {
  OBSOLETE("Quadrature impossible in zero dimensions");
}

template <int dim> Quadrature<dim>::Quadrature(
const Quadrature<dim-1> &q1,const Quadrature<1> &q2) : 
number_points(q1.numberPoints()*q2.numberPoints()) {
  quadrature_points.allocate(number_points);
  quadrature_weights.allocate(number_points);
  int k=0;
  for (int i=0;i<q1.numberPoints();++i) {
    for (int j=0;j<q2.numberPoints();++j) {
      for (int d=0;d<dim-1;++d) {
        quadrature_points[k][d]=q1.point(i)[d];
      }
      quadrature_points[k][dim-1]=q2.point(j)[0];
      quadrature_weights[k]=q1.weight(i)*q2.weight(j);
      ++k;
    };
  }
#ifdef DEBUG
  REAL sum=0.;
  for (int i=0;i<number_points;++i) sum+=quadrature_weights[i];
  CHECK_TEST(abs(sum-1.)<DBL_EPSILON*static_cast<REAL>(number_points));
#endif
}

template <int dim> void Quadrature<dim>::printOn(ostream &os) const {
  os << "Quadrature<" << dim << ">: number_points = " << number_points
     << endl;
  os << "\tquadrature_points:" << endl;
  for (int i=0;i<quadrature_points.getNumber();i++) {
    os << "\t\tquadrature_points[" << i << "] = " 
       << quadrature_points[i] << endl;
  }
  os << "\tquadrature_weights:" << endl;
  for (int i=0;i<quadrature_weights.getNumber();i++) {
    os << "\t\tquadrature_weights[" << i << "] = " 
       << quadrature_weights[i] << endl;
  }
}

template class Quadrature<0>;
template class Quadrature<1>;
#if (SPACEDIM>1)
template class Quadrature<2>;
#endif
#if (SPACEDIM==3)
template class Quadrature<3>;
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int dimensionalRoot(int dim,int n) {
  switch (dim) {
    case 1:
      return n;
    case 2: {
      int i=static_cast<int>(nearbyint(sqrt(static_cast<REAL>(n))));
      CHECK_SAME(i*i,n);
      return i;
    }
    case 3: {
      int i=static_cast<int>(nearbyint(pow(static_cast<REAL>(n),1./3.)));
      CHECK_SAME(i*i*i,n);
      return i;
    }
    default:
      OBSOLETE("case not programmed");
  }
}
#if (SPACEDIM>=2)
int triangleRoot(int n) {
  int i=static_cast<int>(
        nearbyint((sqrt(1.+8.*static_cast<REAL>(n))-1.)*0.5));
  CHECK_SAME(i*(i+1),2*n);
  return i;
}
#endif
#if (SPACEDIM==3)
int tetrahedronRoot(int n) {
  REAL n3=static_cast<REAL>(3*n);
  REAL s=sqrt(n3*n3-1./27.);
  REAL A=pow(n3+s,1./3.);
  REAL B=pow(n3-s,1./3.);
  REAL a=A+B-1.;
  int i=static_cast<int>(nearbyint(a));
  CHECK_SAME(i*(i+1)*(i+2),6*n);
  return i;
}
#endif

extern "C" {
  void F77NAME(dsterf)(const int &n,double *d,double *e,int &info);
}
template <> GaussianQuadrature<1>::GaussianQuadrature(int n) :
Quadrature<1>(n),np_1d(n) {
  CHECK_POSITIVE(n);
  double *d=OPERATOR_NEW_BRACKET(double,n);
  for (int i=0;i<n;i++) d[i]=0.;
  double *e=OPERATOR_NEW_BRACKET(double,n-1);
  for (int i=0;i<n-1;i++) {
    e[i]=static_cast<double>(i+1)
        /sqrt(static_cast<double>((2*i+1)*(2*i+3)));
  }
  int info=0;
  F77NAME(dsterf)(n,d,e,info);
  CHECK_TEST(info==0);
  LegendrePolynomial lp;
  for (int i=0;i<n;i++) {
    REAL xi=(d[i]+1.)*0.5;
    this->quadrature_points[i]=Point<1>(xi);
    NumPtr<REAL> slopes(2);
    lp.derivatives(xi,n,slopes);
    this->quadrature_weights[i]=4./((1.-d[i]*d[i])*slopes[1]*slopes[1]);
  }
  delete [] e;
  delete [] d;

}

template <int dim> void GaussianQuadrature<dim>::printOn(ostream &os) 
const {
  os << "GaussianQuadrature<" << dim << ">: np_1d = " << np_1d << endl;
  Quadrature<dim>::printOn(cout);
}

template class GaussianQuadrature<1>;
#if (SPACEDIM>1)
template <int dim> GaussianQuadrature<dim>::GaussianQuadrature(
int n) : np_1d(n),
  Quadrature<dim>(GaussianQuadrature<dim-1>(n),GaussianQuadrature<1>(n))
  {}
template class GaussianQuadrature<2>;
#endif
#if (SPACEDIM==3)
template class GaussianQuadrature<3>;
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
void F77NAME(dgesv)(const int &n,const int &nrhs,double *a,
  const int &lda,int *ipiv,double *b,const int &ldb,int &info);
}

template <> LobattoQuadrature<1>::LobattoQuadrature(int n) :
Quadrature<1>(n),np_1d(n) {
  CHECK_TEST(n>=2);
  this->quadrature_points[0]=0.;
  this->quadrature_points[n-1]=1.;
  REAL w=this->quadrature_weights[0]=this->quadrature_weights[n-1]
    =1./static_cast<REAL>(n*(n-1));
  if (n>2) {
    double *d=OPERATOR_NEW_BRACKET(double,n-2);
    for (int i=0;i<n-2;i++) d[i]=0.;
    double *e=OPERATOR_NEW_BRACKET(double,n-3);
    for (int i=0;i<n-3;i++) {
      e[i]=sqrt(static_cast<double>((i+1)*(i+3))
               /static_cast<double>((2*i+3)*(2*i+5)));
    }
    int info=0;
    F77NAME(dsterf)(n-2,d,e,info);
    CHECK_TEST(info==0);
    NumPtr<REAL> val(n);
    LegendrePolynomial lp;
    for (int i=1;i<n-1;i++) {
      REAL xi=(d[i-1]+1.)*0.5;
      this->quadrature_points[i]=Point<1>(xi);
      lp.values(xi,val);
      this->quadrature_weights[i]=w/pow(val[n-1],2);
    }
    delete [] e;
    delete [] d;
  }
}

template <int dim> void LobattoQuadrature<dim>::printOn(ostream &os) 
const {
  os << "LobattoQuadrature<" << dim << ">: np_1d = " << np_1d << endl;
  Quadrature<dim>::printOn(cout);
}

template class LobattoQuadrature<1>;
#if (SPACEDIM>1)
template <int dim> LobattoQuadrature<dim>::LobattoQuadrature(
int n) : np_1d(n),
  Quadrature<dim>(LobattoQuadrature<dim-1>(n),LobattoQuadrature<1>(n))
  {}
template class LobattoQuadrature<2>;
#endif
#if (SPACEDIM==3)
template class LobattoQuadrature<3>;
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//see Ralson & Rabinowitz, A First Course in Numerical Analysis, p 118ff
//note that quadrature points are zeros of C^0 Lagrange polynomials
//for n>8 the rules have negative weights
template <> NewtonCotesQuadrature<1>::NewtonCotesQuadrature(int n) : 
Quadrature<1>(n),np_1d(n) {
  CHECK_BOUNDS(n,2,9);
  REAL h=1./static_cast<REAL>(n-1);
  REAL x=0.;
  for (int i=0;i<n;i++,x+=h) {
    this->quadrature_points[i]=x;
  }
  switch (n) {
    case 1:
    case 2: // trapezoidal rule
      this->quadrature_weights[1]=this->quadrature_weights[0]=0.5;
      break;
    case 3: // simpson's rule
      this->quadrature_weights[2]=this->quadrature_weights[0]=1./6.;
      this->quadrature_weights[1]=2./3.;
      break;
    case 4:
      this->quadrature_weights[3]=this->quadrature_weights[0]=0.125;
      this->quadrature_weights[2]=this->quadrature_weights[1]=0.375;
      break;
    case 5:
      this->quadrature_weights[4]=this->quadrature_weights[0]=7./90.;
      this->quadrature_weights[3]=this->quadrature_weights[1]=16./45.;
      this->quadrature_weights[2]=2./15.;
      break;
    case 6:
      this->quadrature_weights[5]=this->quadrature_weights[0]=19./288.;
      this->quadrature_weights[4]=this->quadrature_weights[1]=25./96.;
      this->quadrature_weights[3]=this->quadrature_weights[2]=25./144.;
      break;
    case 7:
      this->quadrature_weights[6]=this->quadrature_weights[0]=41./840.;
      this->quadrature_weights[5]=this->quadrature_weights[1]=216./840.;
      this->quadrature_weights[4]=this->quadrature_weights[2]=27./840.;
      this->quadrature_weights[3]=272./840.;
      break;
    case 8:
    default:
      this->quadrature_weights[7]=this->quadrature_weights[0]
                                 =751./17280.;
      this->quadrature_weights[6]=this->quadrature_weights[1]
                                 =3577./17280.;
      this->quadrature_weights[5]=this->quadrature_weights[2]
                                 =1323./17280.;
      this->quadrature_weights[4]=this->quadrature_weights[3]
                                 =2989./17280.;
      break;
  }
}

template <int dim> void NewtonCotesQuadrature<dim>::printOn(ostream &os)
const {
  os << "NewtonCotesQuadrature<" << dim << ">: np_1d = " << np_1d 
     << endl;
  Quadrature<dim>::printOn(cout);
}

template class NewtonCotesQuadrature<1>;
#if (SPACEDIM>1)
template <int dim> NewtonCotesQuadrature<dim>::NewtonCotesQuadrature(
int n) : Quadrature<dim>(NewtonCotesQuadrature<dim-1>(n),
  NewtonCotesQuadrature<1>(n)),np_1d(n) {}
template class NewtonCotesQuadrature<2>;
#endif
#if (SPACEDIM==3)
template class NewtonCotesQuadrature<3>;
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <> ClenshawCurtisQuadrature<1>::ClenshawCurtisQuadrature(int n)
: Quadrature<1>(n),np_1d(n) {
  CHECK_TEST(n>=2);
  int N=n-1;
  REAL dtheta=0.5*M_PI/static_cast<REAL>(N);
  this->quadrature_points[0]=0.;
  REAL theta=dtheta;
  for (int i=1;i<N;i++,theta+=dtheta) {
    this->quadrature_points[i]=pow(sin(theta),2);
  }
  this->quadrature_points[N]=1.;

  dtheta *= 4.;
  int M=N/2;
  REAL h=1./static_cast<REAL>(N);
  for (int k=0;k<=M;k++) {
    REAL sum=0.;
    for (int j=1;j<M;j++) {
      sum+=cos(static_cast<REAL>(j*k)*dtheta)/static_cast<REAL>(1-4*j*j);
    }
    sum=2.*sum+1.;
    if (M>0) {
      REAL term=cos(M*k*dtheta)/static_cast<REAL>(1.-4*M*M);
      sum+=(N==2*M ? term : 2.*term);
    }
    sum*=h;
    if (k==0) sum*=0.5;
    this->quadrature_weights[N-k]=this->quadrature_weights[k]=sum;
  }
#ifdef DEBUG
//printOn(cout);
#endif
}

template <int dim> void ClenshawCurtisQuadrature<dim>::printOn(
ostream &os) const {
  os << "ClenshawCurtisQuadrature<" << dim << ">: np_1d = " << np_1d 
     << endl;
  Quadrature<dim>::printOn(cout);
}

template class ClenshawCurtisQuadrature<1>;
#if (SPACEDIM>1)
template <int dim>
ClenshawCurtisQuadrature<dim>::ClenshawCurtisQuadrature(int n) :
Quadrature<dim>(ClenshawCurtisQuadrature<dim-1>(n),
  ClenshawCurtisQuadrature<1>(n)),np_1d(n) {}
template class ClenshawCurtisQuadrature<2>;
#endif
#if (SPACEDIM==3)
template class ClenshawCurtisQuadrature<3>;
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
ProductQuadrature::ProductQuadrature(const Quadrature<2> &q2d,
const Quadrature<1> &q1d) {
  base_order=min(q2d.order(),q1d.order());
  number_points=q2d.numberPoints()*q1d.numberPoints();
  quadrature_points.allocate(number_points);
  quadrature_weights.allocate(number_points);
  int ij=0;
  for (int j=0;j<q1d.numberPoints();j++) {
    Point<1> p1d=q1d.point(j);
    REAL w1d=q1d.weight(j);
    for (int i=0;i<q2d.numberPoints();i++,ij++) {
      Point<2> p2d=q2d.point(i);
      quadrature_points[ij]=Point<3>(p2d[0],p2d[1],p1d[0]);
      quadrature_weights[ij]=q2d.weight(i)*w1d;
    }
  }
}
void ProductQuadrature::printOn(ostream &os) const {
  os << "ProductQuadrature: base_order = " << base_order << endl;
  Quadrature<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void set0(NumPtr<Point<2> > &points,NumPtr<REAL> &weights,double w) {
  double xi=1./3.;
  points[0]=Point<2>(xi,xi);
  weights[0]=0.5*w;
  if (!points.ownsData()) ++points;
  if (!weights.ownsData()) ++weights;
}
#if (SPACEDIM==3)
void set0(NumPtr<Point<3> > &points,NumPtr<REAL> &weights,double w) {
  double xi=0.25;
  points[0]=Point<3>(xi,xi,xi);
  weights[0]=w/6.;
  if (!points.ownsData()) ++points;
  if (!weights.ownsData()) ++weights;
}
#endif
void set1(NumPtr<Point<2> > &points,NumPtr<REAL> &weights,
double xi_repeated,double w) {
  double xi_other=1.-2.*xi_repeated;
  points[0]=Point<2>(xi_repeated,xi_repeated);
  points[1]=Point<2>(xi_repeated,xi_other);
  points[2]=Point<2>(xi_other,xi_repeated);
  w*=0.5;
  for (int i=0;i<3;i++) weights[i]=w;
  if (!points.ownsData()) points+=3;
  if (!weights.ownsData()) weights+=3;
}
#if (SPACEDIM==3)
void set1(NumPtr<Point<3> > &points,NumPtr<REAL> &weights,
double xi_repeated,double w) {
  double xi_other=1.-3.*xi_repeated;
  points[0]=Point<3>(xi_repeated,xi_repeated,xi_repeated);
  points[1]=Point<3>(xi_repeated,xi_repeated,xi_other);
  points[2]=Point<3>(xi_repeated,xi_other,xi_repeated);
  points[3]=Point<3>(xi_other,xi_repeated,xi_repeated);
  w/=6.;
  for (int i=0;i<4;i++) weights[i]=w;
  if (!points.ownsData()) points+=4;
  if (!weights.ownsData()) weights+=4;
}
void setPairs(NumPtr<Point<3> > &points,NumPtr<REAL> &weights,
double xi_first_pair,double w) {
  double xi_other_pair=0.5-xi_first_pair;
  points[0]=Point<3>(xi_first_pair,xi_first_pair,xi_other_pair);
  points[1]=Point<3>(xi_first_pair,xi_other_pair,xi_first_pair);
  points[2]=Point<3>(xi_first_pair,xi_other_pair,xi_other_pair);
  points[3]=Point<3>(xi_other_pair,xi_first_pair,xi_first_pair);
  points[4]=Point<3>(xi_other_pair,xi_first_pair,xi_other_pair);
  points[5]=Point<3>(xi_other_pair,xi_other_pair,xi_first_pair);
  w/=6.;
  for (int i=0;i<6;i++) weights[i]=w;
  if (!points.ownsData()) points+=6;
  if (!weights.ownsData()) weights+=6;
}
#endif
void set2(NumPtr<Point<2> > &points,NumPtr<REAL> &weights,
double xi_0,double xi_1,double w) {
  double xi_2=1.-xi_0-xi_1;
  points[0]=Point<2>(xi_0,xi_1);
  points[1]=Point<2>(xi_1,xi_2);
  points[2]=Point<2>(xi_2,xi_0);
  points[3]=Point<2>(xi_2,xi_1);
  points[4]=Point<2>(xi_1,xi_0);
  points[5]=Point<2>(xi_0,xi_2);
  w*=0.5;
  for (int i=0;i<6;i++) weights[i]=w;
  if (!points.ownsData()) points+=6;
  if (!weights.ownsData()) weights+=6;
}
#if (SPACEDIM==3)
void set2(NumPtr<Point<3> > &points,NumPtr<REAL> &weights,
double xi_0,double xi_1,double w) {
  double xi_2=(1.-xi_0-xi_1)*0.5;
  points[0]=Point<3>(xi_0,xi_1,xi_2);
  points[1]=Point<3>(xi_0,xi_2,xi_1);
  points[2]=Point<3>(xi_0,xi_2,xi_2);
  points[3]=Point<3>(xi_1,xi_0,xi_2);
  points[4]=Point<3>(xi_1,xi_2,xi_0);
  points[5]=Point<3>(xi_1,xi_2,xi_2);
  points[6]=Point<3>(xi_2,xi_0,xi_1);
  points[7]=Point<3>(xi_2,xi_0,xi_2);
  points[8]=Point<3>(xi_2,xi_1,xi_0);
  points[9]=Point<3>(xi_2,xi_1,xi_2);
  points[10]=Point<3>(xi_2,xi_2,xi_0);
  points[11]=Point<3>(xi_2,xi_2,xi_1);
  w/=6.;
  for (int i=0;i<12;i++) weights[i]=w;
  if (!points.ownsData()) points+=12;
  if (!weights.ownsData()) weights+=12;
}
void set3(NumPtr<Point<3> > &points,NumPtr<REAL> &weights,
double xi_0,double xi_1,double xi_2,double w) {
  double xi_3=1.-xi_0-xi_1-xi_2;
  points[ 0]=Point<3>(xi_0,xi_1,xi_2);
  points[ 1]=Point<3>(xi_0,xi_1,xi_3);
  points[ 2]=Point<3>(xi_0,xi_2,xi_1);
  points[ 3]=Point<3>(xi_0,xi_2,xi_3);
  points[ 4]=Point<3>(xi_0,xi_3,xi_1);
  points[ 5]=Point<3>(xi_0,xi_3,xi_2);
  points[ 6]=Point<3>(xi_1,xi_0,xi_2);
  points[ 7]=Point<3>(xi_1,xi_0,xi_3);
  points[ 8]=Point<3>(xi_1,xi_2,xi_0);
  points[ 9]=Point<3>(xi_1,xi_2,xi_3);
  points[10]=Point<3>(xi_1,xi_3,xi_0);
  points[11]=Point<3>(xi_1,xi_3,xi_2);
  points[12]=Point<3>(xi_2,xi_0,xi_1);
  points[13]=Point<3>(xi_2,xi_0,xi_3);
  points[14]=Point<3>(xi_2,xi_1,xi_0);
  points[15]=Point<3>(xi_2,xi_1,xi_3);
  points[16]=Point<3>(xi_2,xi_3,xi_0);
  points[17]=Point<3>(xi_2,xi_3,xi_1);
  points[18]=Point<3>(xi_3,xi_0,xi_1);
  points[19]=Point<3>(xi_3,xi_0,xi_2);
  points[20]=Point<3>(xi_3,xi_1,xi_0);
  points[21]=Point<3>(xi_3,xi_1,xi_2);
  points[22]=Point<3>(xi_3,xi_2,xi_0);
  points[23]=Point<3>(xi_3,xi_2,xi_1);
  w/=6.;
  for (int i=0;i<24;i++) weights[i]=w;
  if (!points.ownsData()) points+=24;
  if (!weights.ownsData()) weights+=24;
}
#endif
// Peter Silvester,
// Symmetric Quadrature Formulae for Simplexes,
// Mathematics of Computation,
// Volume 24, Number 109, January 1970, pages 95-100.
TriangleNewtonCotesQuadrature::TriangleNewtonCotesQuadrature(int np_1d){
  switch (np_1d) {
    case 1:
    case 2: {
      base_order=2;
      number_points=3;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,0.,1./3.);
      break;
    }
    case 3: {
      base_order=3;
      number_points=3;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,0.5,1./3.);
      break;
    }
    case 4: {
      base_order=3;
      number_points=1+3; // = 4
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.75);
      set1(qp,qw,0.,1./12.);
      break;
    }
    case 5: {
      base_order=4;
      number_points=3*3; // = 9
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,2./45.);
      set1(qp,qw,0.5,1./9.);
      set1(qp,qw,0.25,8./45.);
      break;
    }
    case 6:
    default: {
      base_order=4;
      number_points=1+2*3; // = 7
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.45);
      set1(qp,qw,0.,0.05);
      set1(qp,qw,0.5,2./15.);
      break;
    }
  }
}
void TriangleNewtonCotesQuadrature::printOn(ostream &os) const {
  os << "TriangleNewtonCotesQuadrature: base_order = " << base_order
     << endl;
  Quadrature<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
TriangleFullNewtonCotesQuadrature::TriangleFullNewtonCotesQuadrature(
int np_1d) {
  switch (np_1d) {
    case 2: {
      base_order=2;
      number_points=3;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,0.,1./3.);
      break;
    }
    case 3: {
      base_order=3;
      number_points=2*3; // = 6
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.5,1./3.);
      set1(qp,qw,0.,0.);
      break;
    }
    case 4: {
      base_order=3;
      number_points=1+3+6; // = 10
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.75);
      set1(qp,qw,0.,1./12.);
      set2(qp,qw,0.,1./3.,0.);
      break;
    }
    case 5: {
//
      base_order=4;
      number_points=3*3+6; // = 15
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,2./45.);
      set1(qp,qw,0.5,1./9.);
      set1(qp,qw,0.25,8./45.);
      set2(qp,qw,0.,0.25,0.);
      break;
    }
    case 6: {
      base_order=4;
      number_points=3*3+2*6; // = 21
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,11./1008.);
      set1(qp,qw,0.2,200./1008.);
      set1(qp,qw,0.4,25./1008.);
      set2(qp,qw,0.,0.2,25./1008.);
      set2(qp,qw,0.,0.4,25./1008.);
      break;
    }
    case 7: {
//
      base_order=4;
      number_points=1+3*3+3*6; // = 28
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.45);
      set1(qp,qw,0.,0.05);
      set1(qp,qw,0.5,2./15.);
      set1(qp,qw,1./6.,0.);
      set2(qp,qw,0.,1./6.,0.);
      set2(qp,qw,0.,1./3.,0.);
      set2(qp,qw,1./6.,1./3.,0.);
      break;
    }
    case 8: {
      base_order=5;
      number_points=4*3+4*6; // = 36
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,1336./259200.);
      set1(qp,qw,1./7.,32242./259200.);
      set1(qp,qw,2./7.,3430./259200.);
      set1(qp,qw,3./7.,44590./259200.);
      set2(qp,qw,0.,1./7.,2989./259200.);
      set2(qp,qw,0.,2./7.,3577./259200.);
      set2(qp,qw,0.,3./7.,2695./259200.);
      set2(qp,qw,1./7.,2./7.,-6860./259200.);
      break;
    }
    case 9: 
    default: { 
      base_order=5;
      number_points=1+3*3+3*6; // = 28
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,0.);
      set1(qp,qw,0.125,704./14175.);
      set1(qp,qw,0.25,-1448./14175.);
      set1(qp,qw,0.375,1472./14175.);
      set1(qp,qw,0.5,-1083./14175.);
      set2(qp,qw,0.,0.125,368./14175.);
      set2(qp,qw,0.,0.25,-468/14175.);
      set2(qp,qw,0.,0.375,1136./14175.);
      set2(qp,qw,0.125,0.25,832./14175.);
      set2(qp,qw,0.125,0.375,672./14175.);
      break;
    }
  }
}
void TriangleFullNewtonCotesQuadrature::printOn(ostream &os) const {
  os << "TriangleFullNewtonCotesQuadrature: base_order = " << base_order
     << endl;
  Quadrature<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
TriangleGaussianQuadrature::TriangleGaussianQuadrature(int ord) {
  switch (ord) {
    case 0:
    case 1:
    case 2: // dunavent 1
      base_order=2;
      number_points=1;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set0(quadrature_points,quadrature_weights,1.);
      break;
    case 3: {
//    exact for | alpha | <= 2
      base_order=3;
      number_points=3;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,0.5,1./3.);
      break;
    }
    case 4: 
    case 5: { // dunavant 4
  //  exact for | alpha | <=4
      base_order=5;
      number_points=2*3; // = 6
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
//    b1 is root of 135 z^4 - 240 z^3 + 120 z^2 - 20 z + 1
//    b0=-75 b1^3 + 355/3 b1^2 - 43 b1 +28/9
//    w0=-585/62 b1^3 + 3935/248 b1^2 - 807/124 b1 + 43/62
//    w1-1-3 - w0
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,.091576213509770743460,0.10995174365532186765);
      set1(qp,qw,.44594849091596488635,0.22338158967801146568);
      break;
    }
    case 6: { // wandzura 1
//    exact for | alpha | <= 5
      base_order=6;
      number_points=1+2*3; // = 7
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,.225);
//    b1, c1 solve 21 z^2 - 12 z + 1 = 0
      set1(qp,qw,(6.+sqrt(15.))/21.,(155.+sqrt(15.))/1200);
      set1(qp,qw,(6.-sqrt(15.))/21.,(155.-sqrt(15.))/1200);
      break;
    }
    case 7: // dunavant 7 is unstable
    case 8: 
    case 9:{ // dunavant 8
      number_points=1+3*3+6; // = 16
      base_order=9;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.144315607677787);
      set1(qp,qw,0.459292588292723,0.095091634267285);
      set1(qp,qw,0.170569307751760,0.103217370534718);
      set1(qp,qw,0.050547228317031,0.032458497623198);
      set2(qp,qw,0.008394777409958,0.263112829634638,
        0.027230314174435);
      break;
    }
    case 10: { // dunavant 9
      number_points=1+4*3+6; // = 19
      base_order=10;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.097135796282799);
      set1(qp,qw,0.489682519198738,0.031334700227139);
      set1(qp,qw,0.437089591492937,0.077827541004774);
      set1(qp,qw,0.188203535619033,0.079647738927210);
      set1(qp,qw,0.044729513394453,0.025577675658698);
      set2(qp,qw,0.036838412054736,0.221962989160766,
        0.043283539377289);
      break;
    }
//  dunavant 10 is inefficient
    case 11: { // wandzura 2
      base_order=11;
      number_points=1+4*3+2*6; // = 25
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.8352339980519638E-01);
      set1(qp,qw,0.49786543295447,0.7229850592056743E-02);
      set1(qp,qw,0.42801244972906,0.7449217792098051E-01);
      set1(qp,qw,0.18475641274322,0.7864647340310853E-01);
      set1(qp,qw,0.02048121857168,0.6928323087107504E-02);
      set2(qp,qw,0.03500298989727,0.13657357625603,
        0.2951832033477940E-01);
      set2(qp,qw,0.03754907025844,0.33274360058864,
        0.3957936719606124E-01);
      break;
    }
    case 12:
    case 13: { // dunavant 12
      number_points=5*3+3*6; // = 33
      base_order=13;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.488217389773805,0.025731066440455);
      set1(qp,qw,0.439724392294460,0.043692544538038);
      set1(qp,qw,0.271210385012116,0.062858224217885);
      set1(qp,qw,0.127576145541586,0.034796112930709);
      set1(qp,qw,0.021317350453210,0.006166261051559);
      set2(qp,qw,0.115343494534698,0.275713269685514,
        0.040371557766381);
      set2(qp,qw,0.022838332222257,0.281325580989940,
        0.022356773202303);
      set2(qp,qw,0.025734050548330,0.116251915907597,
        0.017316231108659);
      break;
    }
    case 14: { // dunavant 13
      number_points=1+6*3+3*6; // = 37
      base_order=14;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.052520923400802);
      set1(qp,qw,0.495048184939705,0.011280145209330);
      set1(qp,qw,0.468716635109574,0.031423518362454);
      set1(qp,qw,0.414521336801277,0.047072502504194);
      set1(qp,qw,0.229399572042831,0.047363586536355);
      set1(qp,qw,0.114424495196330,0.031167529045794);
      set1(qp,qw,0.024811391363459,0.007975771465074);
      set2(qp,qw,0.094853828379579,0.268794997058761,
        0.036848402728732);
      set2(qp,qw,0.018100773278807,0.291730066734288,
        0.017401463303822);
      set2(qp,qw,0.022233076674090,0.126357385491669,
        0.015521786839045);
      break;
    }
    case 15: { // dunavant 14
      number_points=6*3+4*6; // = 42
      base_order=15;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.488963910362179,0.021883581369429);
      set1(qp,qw,0.417644719340454,0.032788353544125);
      set1(qp,qw,0.273477528308839,0.051774104507292);
      set1(qp,qw,0.177205532412543,0.042162588736993);
      set1(qp,qw,0.061799883090873,0.014433699669777);
      set1(qp,qw,0.019390961248701,0.004923403602400);
      set2(qp,qw,0.057124757403648,0.172266687821356,
        0.024665753212564);
      set2(qp,qw,0.092916249356972,0.336861459796345,
        0.038571510787061);
      set2(qp,qw,0.014646950055654,0.298372882136258,
        0.014436308113534);
      set2(qp,qw,0.001268330932872,0.118974497696957,
        0.005010228838501);
      break;
    }
//  dunavant 15 is unstable 
    case 16: { // wandzura 3
      base_order=16;
      number_points=6*3+6*6; // = 54
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.45828079636912,0.3266181884880529E-01);
      set1(qp,qw,0.40361046457913,0.2741281803136436E-01);
      set1(qp,qw,0.29319716791303,0.2651003659870330E-01);
      set1(qp,qw,0.14646778694277,0.2921596213648611E-01);
      set1(qp,qw,0.05636286766560,0.1058460806624399E-01);
      set1(qp,qw,0.01657512685837,0.3614643064092035E-02);
      set2(qp,qw,0.00991220330923,0.23953455415479,
        0.8527748101709436E-02);
      set2(qp,qw,0.01580377063023,0.40487880731834,
        0.1391617651669193E-01);
      set2(qp,qw,0.00514360881697,0.09500211311304,
        0.4291932940734835E-02);
      set2(qp,qw,0.04892232575299,0.14975310732227,
        0.1623532928177489E-01);
      set2(qp,qw,0.06876874863252,0.28691961244133,
        0.2560734092126239E-01);
      set2(qp,qw,0.16840441812470,0.28183566809908,
        0.3308819553164567E-01);
      break;
    }
    case 17:
    case 18: { // dunavant 17
      number_points=1+8*3+6*6; // = 61
      base_order=18;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.033437199290803);
      set1(qp,qw,0.497170540556774,0.005093415440507);
      set1(qp,qw,0.482176322624625,0.014670864527638);
      set1(qp,qw,0.450239969020782,0.024350878353672);
      set1(qp,qw,0.400266239377397,0.031107550868969);
      set1(qp,qw,0.252141267970953,0.031257111218620);
      set1(qp,qw,0.162047004658461,0.024815654339665);
      set1(qp,qw,0.075875882260746,0.014056073070557);
      set1(qp,qw,0.015654726967822,0.003194676173779);
      set2(qp,qw,0.010186928826919,0.334319867363658,
        0.008119655318993);
      set2(qp,qw,0.135440871671036,0.292221537796944,
        0.026805742283163);
      set2(qp,qw,0.054423924290583,0.319574885423190,
        0.018459993210822);
      set2(qp,qw,0.012868560833637,0.190704224192292,
        0.008476868534328);
      set2(qp,qw,0.067165782413524,0.180483211648746,
        0.018292796770025);
      set2(qp,qw,0.014663182224828,0.080711313679564,
        0.006665632004165);
      break;
    }
    case 19:
    case 20: { // dunavant 19
      number_points=1+8*3+8*6; // = 73
      base_order=20;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.032906331388919);
      set1(qp,qw,0.489609987073006,0.010330731891272);
      set1(qp,qw,0.454536892697893,0.022387247263016);
      set1(qp,qw,0.401416680649431,0.030266125869468);
      set1(qp,qw,0.255551654403098,0.030490967802198);
      set1(qp,qw,0.177077942152130,0.024159212741641);
      set1(qp,qw,0.110061053227952,0.016050803586801);
      set1(qp,qw,0.055528624251840,0.008084580261784);
      set1(qp,qw,0.012621863777229,0.002079362027485);
      set2(qp,qw,0.003611417848412,0.395754787356943,
        0.003884876904981);
      set2(qp,qw,0.134466754530780,0.307929983880436,
        0.025574160612022);
      set2(qp,qw,0.014446025776115,0.264566948406520,
        0.008880903573338);
      set2(qp,qw,0.046933578838178,0.358539352205951,
        0.016124546761731);
      set2(qp,qw,0.002861120350567,0.157807405968595,
        0.002491941817491);
      set2(qp,qw,0.223861424097916,0.075050596975911,
        0.018242840118951);
      set2(qp,qw,0.034647074816760,0.142421601113383,
        0.010258563736199);
      set2(qp,qw,0.010161119296278,0.065494628082938,
        0.003799928855302);
      break;
    }
//  dunavant 20 is unstable
    case 21: { // wandzura 4
      base_order=21;
      number_points=1+8*3+10*6; // 85
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.2761042699769952E-01);
      set1(qp,qw,0.49924967533779,0.1779029547326740E-02);
      set1(qp,qw,0.45293012403052,0.2011239811396117E-01);
      set1(qp,qw,0.39776393795524,0.2681784725933157E-01);
      set1(qp,qw,0.26450020253279,0.2452313380150201E-01);
      set1(qp,qw,0.21101896409208,0.1639457841069539E-01);
      set1(qp,qw,0.10773560717127,0.1479590739864960E-01);
      set1(qp,qw,0.03906908783780,0.4579282277704251E-02);
      set1(qp,qw,0.01117437972933,0.1651826515576217E-02);
      set2(qp,qw,0.00534961818734,0.06354966590835,
        0.2349170908575584E-02);
      set2(qp,qw,0.00795481706620,0.15710691894071,
        0.4465925754181793E-02);
      set2(qp,qw,0.01042239828126,0.39564211436437,
        0.6099566807907972E-02);
      set2(qp,qw,0.01096441479612,0.27316757071291,
        0.6891081327188203E-02);
      set2(qp,qw,0.03856671208546,0.10178538248502,
        0.7997475072478163E-02);
      set2(qp,qw,0.03558050781722,0.44665854917641,
        0.7386134285336024E-02);
      set2(qp,qw,0.04967081636276,0.19901079414950,
        0.1279933187864826E-01);
      set2(qp,qw,0.05851972508433,0.32426118369228,
        0.1725807117569655E-01);
      set2(qp,qw,0.12149778700439,0.20853136321013,
        0.1867294590293547E-01);
      set2(qp,qw,0.14071084494394,0.32317056653626,
        0.2281822405839526E-01);
      break;
    }
    case 22:
    case 23:
    case 24:
    case 25:
    case 26: { // wandzura 5
      base_order=26;
      number_points=10*3+16*6; // 126
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.48602675846341,0.8005581880020417E-02);
      set1(qp,qw,0.43441069933617,0.1594707683239050E-01);
      set1(qp,qw,0.38988913524396,0.1310914123079553E-01);
      set1(qp,qw,0.29844323401980,0.1958300096563562E-01);
      set1(qp,qw,0.23404417233737,0.1647088544153727E-01);
      set1(qp,qw,0.15146833460902,0.8547279074092100E-02);
      set1(qp,qw,0.11273389354599,0.8161885857226492E-02);
      set1(qp,qw,0.07771569209153,0.6121146539983779E-02);
      set1(qp,qw,0.03489309361430,0.2908498264936665E-02);
      set1(qp,qw,0.00725818462093,0.6922752456619963E-03);
      set2(qp,qw,0.00129235270444,0.22721445215336,
        0.1248289199277397E-02);
      set2(qp,qw,0.00539970127212,0.43501055485357,
        0.3404752908803022E-02);
      set2(qp,qw,0.00638400303398,0.32030959927220,
        0.3359654326064051E-02);
      set2(qp,qw,0.00502821150199,0.09175032228001,
        0.1716156539496754E-02);
      set2(qp,qw,0.00682675862178,0.03801083585872,
        0.1480856316715606E-02);
      set2(qp,qw,0.01001619963993,0.15742521848531,
        0.3511312610728685E-02);
      set2(qp,qw,0.02575781317339,0.23988965977853,
        0.7393550149706484E-02);
      set2(qp,qw,0.03022789811992,0.36194311812606,
        0.7983087477376558E-02);
      set2(qp,qw,0.03050499010716,0.08355196095483,
        0.4355962613158041E-02);
      set2(qp,qw,0.04595654736257,0.14844322073242,
        0.7365056701417832E-02);
      set2(qp,qw,0.06744280054028,0.28373970872753,
        0.1096357284641955E-01);
      set2(qp,qw,0.07004509141591,0.40689937511879,
        0.1174996174354112E-01);
      set2(qp,qw,0.08391152464012,0.19411398702489,
        0.1001560071379857E-01);
      set2(qp,qw,0.12037553567715,0.32413434700070,
        0.1330964078762868E-01);
      set2(qp,qw,0.14806689915737,0.22927748355598,
        0.1415444650522614E-01);
      set2(qp,qw,0.19177186586733,0.32561812259598,
        0.1488137956116801E-01);
      break;
    }
    case 27:
    case 28:
    case 29:
    case 30:
    case 31: 
    default: { // wandzura 6
      base_order=31;
      number_points=1+12*3+23*6; // 175
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.1557996020289920E-01);
      set1(qp,qw,0.49633494178362,0.3177233700534134E-02);
      set1(qp,qw,0.45850216209852,0.1048342663573077E-01);
      set1(qp,qw,0.42450952193729,0.1320945957774363E-01);
      set1(qp,qw,0.38204707005392,0.1497500696627150E-01);
      set1(qp,qw,0.28098784579608,0.1498790444338419E-01);
      set1(qp,qw,0.22734897585403,0.1333886474102166E-01);
      set1(qp,qw,0.17455911150873,0.1088917111390201E-01);
      set1(qp,qw,0.12325842720144,0.8189440660893461E-02);
      set1(qp,qw,0.08008422889220,0.5575387588607785E-02);
      set1(qp,qw,0.04777446740790,0.3191216473411976E-02);
      set1(qp,qw,0.02172051468014,0.1296715144327045E-02);
      set1(qp,qw,0.00476467761544,0.2982628261349172E-03);
      set2(qp,qw,0.00092537119335,0.41529527091331,
        0.9989056850788964E-03);
      set2(qp,qw,0.00138592585556,0.06118990978535,
        0.4628508491732533E-03);
      set2(qp,qw,0.00368241545591,0.16490869013691,
        0.1234451336382413E-02);
      set2(qp,qw,0.00390322342416,0.02503506223200,
        0.5707198522432062E-03);
      set2(qp,qw,0.00323324815501,0.30606446515110,
        0.1126946125877624E-02);
      set2(qp,qw,0.00646743211224,0.10707328373022,
        0.1747866949407337E-02);
      set2(qp,qw,0.00324747549133,0.22995754934558,
        0.1182818815031657E-02);
      set2(qp,qw,0.00867509080675,0.33703663330578,
        0.1990839294675034E-02);
      set2(qp,qw,0.01559702646731,0.05625657618206,
        0.1900412795035980E-02);
      set2(qp,qw,0.01797672125369,0.40245137521240,
        0.4498365808817451E-02);
      set2(qp,qw,0.01712424535389,0.24365470201083,
        0.3478719460274719E-02);
      set2(qp,qw,0.02288340534658,0.16538958561453,
        0.4102399036723953E-02);
      set2(qp,qw,0.03273759728777,0.09930187449585,
        0.4021761549744162E-02);
      set2(qp,qw,0.03382101234234,0.30847833306905,
        0.6033164660795066E-02);
      set2(qp,qw,0.03554761446002,0.46066831859211,
        0.3946290302129598E-02);
      set2(qp,qw,0.05053979030687,0.21881529945393,
        0.6644044537680268E-02);
      set2(qp,qw,0.05701471491573,0.37920955156027,
        0.8254305856078458E-02);
      set2(qp,qw,0.06415280642120,0.14296081941819,
        0.6496056633406411E-02);
      set2(qp,qw,0.08050114828763,0.28373128210592,
        0.9252778144146602E-02);
      set2(qp,qw,0.10436706813453,0.19673744100444,
        0.9164920726294280E-02);
      set2(qp,qw,0.11384489442875,0.35588914121166,
        0.1156952462809767E-01);
      set2(qp,qw,0.14536348771552,0.25981868535191,
        0.1176111646760917E-01);
      set2(qp,qw,0.18994565282198,0.32192318123130,
        0.1382470218216540E-01);
      break;
    }
  }
}
void TriangleGaussianQuadrature::printOn(ostream &os) const {
  os << "TriangleGaussianQuadrature: base_order = " << base_order
     << endl;
  Quadrature<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
// mapping is
//   | x_0 | = | xi_0 ( 1 - xi_1 / 2 ) |
//   | x_1 | = | xi_1 ( 1 - xi_0 / 2 ) |
// so
//     d x  = | 1 - xi_1 / 2 ,   - xi_0 / 2 |
//   / d xi = |   - xi_1 / 2 , 1 - xi_0 / 2 |
// and
//   det(d x / d xi ) = 1 - ( xi_0 + xi-1 ) / 2
//mapping determinant is linear, 
//  so Gaussian quad on rectangle is applied to integrand times a linear
//  so order is one less than Gaussian quadrature on a rectangle
TriangleProductGaussianQuadrature::TriangleProductGaussianQuadrature(
int n) {
  GaussianQuadrature<2> gq(n); // ==> order=2*n
  number_points=gq.numberPoints();
  base_order=gq.order()-1; // = 2*n-1
  quadrature_points.allocate(number_points);
  quadrature_weights.allocate(number_points);
  for (int p=0;p<number_points;p++) {
    Point<2> &qp=quadrature_points[p];
    const Point<2> &gq_qp=gq.points()[p];
    qp[0]=gq_qp[0]*(1.-0.5*gq_qp[1]);
    qp[1]=gq_qp[1]*(1.-0.5*gq_qp[0]);
    quadrature_weights[p]=
      gq.weights()[p]*(1.-0.5*(gq_qp[0]+gq_qp[1]));
  }
}
void TriangleProductGaussianQuadrature::printOn(ostream &os) const {
  os << "TriangleProductGaussianQuadrature: base_order = " 
     << base_order << endl;
  Quadrature<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
TriangleLobattoQuadrature::TriangleLobattoQuadrature(int ord) {
  switch (ord) {
    case 0:
    case 1:
    case 2:
//    exact for | alpha | <= 1
      base_order=2;
      number_points=3;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,0.,1./3.);
      break;
    case 3: {
//    exact for | alpha | <= 2
      base_order=3;
      number_points=1+3;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.75);
      set1(qp,qw,0.,1./12.);
      break;
    }
    case 4: { // fekete 1
      number_points=1+3+6; // = 10
      base_order=4;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.45);
      set1(qp,qw,0.,1./60.);
      set2(qp,qw,0.,(5.-sqrt(5.))/10.,1./12.);
      break;
    }
    case 5:
    case 6:
    case 7: { // fekete 2
      number_points=1+3*3+3*6; // = 28
      base_order=7;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.2178563571);
      set1(qp,qw,0.1063354684,0.5*0.1104193374);
      set1(qp,qw,0.5000000000,0.5*0.0358939762);
      set1(qp,qw,0.0000000000,0.5*0.0004021278);
      set2(qp,qw,0.1171809171,0.3162697959,0.5*0.1771348660);
      set2(qp,qw,0.,0.2655651402,0.5*0.0272344079);
      set2(qp,qw,0.,0.0848854223,0.5*0.0192969460);
      break;
    }
    case 8:
    case 9:
    case 10: { // fekete 3
      number_points=1+4*3+7*6; // = 55
      base_order=10;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.1096011288);
      set1(qp,qw,0.1704318201,0.5*0.0767491008);
      set1(qp,qw,0.4699587644,0.5*0.0646677819);
      set1(qp,qw,0.0489345696,0.5*0.0276211659);
      set1(qp,qw,0.,0.5*0.0013925011);
      set2(qp,qw,0.1784337588,0.3252434900,0.5*0.0933486453);
      set2(qp,qw,0.0588564879,0.3010242110,0.5*0.0619010169);
      set2(qp,qw,0.0551758079,0.1543901944,0.5*0.0437466450);
      set2(qp,qw,0.,0.4173602935,0.5*0.0114553907);
      set2(qp,qw,0.,0.2610371960,0.5*0.0093115568);
      set2(qp,qw,0.,0.1306129092,0.5*0.0078421987);
      set2(qp,qw,0.,0.0402330070,0.5*0.0022457501);
      break;
    }
    case 11:
    case 12:
    case 13: { // fekete 3
      number_points=1+6*3+12*6; // = 91
      base_order=13;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.0626245179);
      set1(qp,qw,0.4005558262,0.5*0.0571359417);
      set1(qp,qw,0.2618405201,0.5*0.0545982307);
      set1(qp,qw,0.0807386775,0.5*0.0172630326);
      set1(qp,qw,0.0336975736,0.5*0.0142519606);
      set1(qp,qw,0.5,0.5*0.0030868485);
      set1(qp,qw,0.,0.5*0.0004270742);
      set2(qp,qw,0.1089969290,0.3837518758,0.5*0.0455876390);
      set2(qp,qw,0.1590834479,0.2454317980,0.5*0.0496701966);
      set2(qp,qw,0.0887037176,0.1697134458,0.5*0.0387998322);
      set2(qp,qw,0.0302317829,0.4071849276,0.5*0.0335323983);
      set2(qp,qw,0.0748751152,0.2874821712,0.5*0.0268431561);
      set2(qp,qw,0.0250122615,0.2489279690,0.5*0.0237377452);
      set2(qp,qw,0.0262645218,0.1206826354,0.5*0.0177255972);
      set2(qp,qw,0.,0.3753565349,0.5*0.0043097313);
      set2(qp,qw,0.,0.2585450895,0.5*0.0028258057);
      set2(qp,qw,0.,0.1569057655,0.5*0.0030994935);
      set2(qp,qw,0.,0.0768262177,0.5*0.0023829062);
      set2(qp,qw,0.,0.0233450767,0.5*0.0009998683);
      break;
    }
    case 14: 
    case 15: 
    case 16: { // fekete 6
      number_points=1+7*3+19*6; // = 136
      base_order=16;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.0459710878);
      set1(qp,qw,0.4206960976,0.5*0.0384470625);
      set1(qp,qw,0.2260541354,0.5*0.0386013566);
      set1(qp,qw,0.1186657611,0.5*0.0224308157);
      set1(qp,qw,0.4761452137,0.5*0.0243531004);
      set1(qp,qw,0.0531173538,0.5*0.0094392654);
      set1(qp,qw,0.0219495841,0.5*0.0061105652);
      set1(qp,qw,0.,0.0001283162);
      set2(qp,qw,0.2379370518,0.3270403780,0.5*0.0346650571);
      set2(qp,qw,0.1585345951,0.3013819154,0.5*0.0305412307);
      set2(qp,qw,0.0972525649,0.3853507643,0.5*0.0262101254);
      set2(qp,qw,0.0875150140,0.2749910734,0.5*0.0265367617);
      set2(qp,qw,0.1339547708,0.1975591066,0.5*0.0269859772);
      set2(qp,qw,0.0475622627,0.3524012205,0.5*0.0172635676);
      set2(qp,qw,0.0596194677,0.1978887556,0.5*0.0188795851);
      set2(qp,qw,0.0534939782,0.1162464503,0.5*0.0158224870);
      set2(qp,qw,0.0157189888,0.4176001732,0.5*0.0127170850);
      set2(qp,qw,0.0196887324,0.2844332752,0.5*0.0164489660);
      set2(qp,qw,0.0180698489,0.1759511193,0.5*0.0120018620);
      set2(qp,qw,0.0171941515,0.0816639421,0.5*0.0072268907);
      set2(qp,qw,0.,0.4493368632,0.5*0.0023599161);
      set2(qp,qw,0.,0.3500847655,0.5*0.0017624674);
      set2(qp,qw,0.,0.2569702891,0.5*0.0018648017);
      set2(qp,qw,0.,0.1738056486,0.5*0.0012975716);
      set2(qp,qw,0.,0.1039958541,0.5*0.0018506035);
      set2(qp,qw,0.,0.0503997335,0.5*0.0009919379);
      set2(qp,qw,0.,0.0152159769,0.5*0.0004893506);
      break;
    }
    case 17: 
    case 18: 
    case 19: 
    default: { // fekete 7
      number_points=1+9*3+27*6; // = 190
      base_order=19;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.0326079297);
      set1(qp,qw,0.4099034502,0.5*0.0288093886);
      set1(qp,qw,0.2438647767,0.5*0.0279490452);
      set1(qp,qw,0.1512564554,0.5*0.0174438045);
      set1(qp,qw,0.4594655253,0.5*0.0203594338);
      set1(qp,qw,0.0832757649,0.5*0.0113349170);
      set1(qp,qw,0.0369065587,0.5*0.0046614185);
      set1(qp,qw,0.0149574850,0.5*0.0030346239);
      set1(qp,qw,0.5,0.5*0.0012508731);
      set1(qp,qw,0.,0.5*0.0000782945);
      set2(qp,qw,0.2515553103,0.3292984162,0.5*0.0255331366);
      set2(qp,qw,0.1821465920,0.3095465041,0.5*0.0235716330);
      set2(qp,qw,0.1246901255,0.3789288931,0.5*0.0206304700);
      set2(qp,qw,0.1179441386,0.2868915642,0.5*0.0204028340);
      set2(qp,qw,0.1639418454,0.2204868669,0.5*0.0215105697);
      set2(qp,qw,0.0742549663,0.3532533654,0.5*0.0183482070);
      set2(qp,qw,0.0937816771,0.2191980979,0.5*0.0174161032);
      set2(qp,qw,0.0890951387,0.1446273457,0.5*0.0155972434);
      set2(qp,qw,0.0409065243,0.4360543636,0.5*0.0119269616);
      set2(qp,qw,0.0488675890,0.2795984854,0.5*0.0147074804);
      set2(qp,qw,0.0460342127,0.2034211147,0.5*0.0116182830);
      set2(qp,qw,0.0420687187,0.1359040280,0.5*0.0087639138);
      set2(qp,qw,0.0116377940,0.4336892286,0.5*0.0098563528);
      set2(qp,qw,0.0299062187,0.3585587824,0.5*0.0096342355);
      set2(qp,qw,0.0132313129,0.2968103667,0.5*0.0086477936);
      set2(qp,qw,0.0136098469,0.2050279257,0.5*0.0083868302);
      set2(qp,qw,0.0124869684,0.1232146223,0.5*0.0062576643);
      set2(qp,qw,0.0365197797,0.0805854893,0.5*0.0077839825);
      set2(qp,qw,0.0118637765,0.0554881302,0.5*0.0031415239);
      set2(qp,qw,0.,0.4154069883,0.5*0.0006513246);
      set2(qp,qw,0.,0.3332475761,0.5*0.0021137942);
      set2(qp,qw,0.,0.2558853572,0.5*0.0004393452);
      set2(qp,qw,0.,0.1855459314,0.5*0.0013662119);
      set2(qp,qw,0.,0.1242528987,0.5*0.0003331251);
      set2(qp,qw,0.,0.0737697111,0.5*0.0011613225);
      set2(qp,qw,0.,0.0355492359,0.5*0.0004342867);
      set2(qp,qw,0.,0.0106941169,0.5*0.0002031499);
      break;
    }
  }
}
void TriangleLobattoQuadrature::printOn(ostream &os) const {
  os << "TriangleLobattoQuadrature: base_order = " << base_order
     << endl;
  Quadrature<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
// Stephen Wandzura, Hong Xiao,
// Symmetric Quadrature Rules on a Triangle,
// Computers and Mathematics with Applications,
// Volume 45, Number 12, June 2003, pages 1829-1840.

TriangleWandzuraQuadrature::TriangleWandzuraQuadrature(int n) {
  switch (n) {
    case 1: {
      base_order=6;
      number_points=1+2*3; // = 7
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.225);
      set1(qp,qw,0.47014206410512,0.1323941527885062);
      set1(qp,qw,0.10128650732346,0.1259391805448271);
      break;
    }
    case 2: {
      base_order=11;
      number_points=1+4*3+2*6; // = 25
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.8352339980519638E-01);
      set1(qp,qw,0.49786543295447,0.7229850592056743E-02);
      set1(qp,qw,0.42801244972906,0.7449217792098051E-01);
      set1(qp,qw,0.18475641274322,0.7864647340310853E-01);
      set1(qp,qw,0.02048121857168,0.6928323087107504E-02);
      set2(qp,qw,0.03500298989727,0.13657357625603,
        0.2951832033477940E-01);
      set2(qp,qw,0.03754907025844,0.33274360058864,
        0.3957936719606124E-01);
      break;
    }
    case 3: {
      base_order=16;
      number_points=6*3+6*6; // = 54
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.45828079636912,0.3266181884880529E-01);
      set1(qp,qw,0.40361046457913,0.2741281803136436E-01);
      set1(qp,qw,0.29319716791303,0.2651003659870330E-01);
      set1(qp,qw,0.14646778694277,0.2921596213648611E-01);
      set1(qp,qw,0.05636286766560,0.1058460806624399E-01);
      set1(qp,qw,0.01657512685837,0.3614643064092035E-02);
      set2(qp,qw,0.00991220330923,0.23953455415479,
        0.8527748101709436E-02);
      set2(qp,qw,0.01580377063023,0.40487880731834,
        0.1391617651669193E-01);
      set2(qp,qw,0.00514360881697,0.09500211311304,
        0.4291932940734835E-02);
      set2(qp,qw,0.04892232575299,0.14975310732227,
        0.1623532928177489E-01);
      set2(qp,qw,0.06876874863252,0.28691961244133,
        0.2560734092126239E-01);
      set2(qp,qw,0.16840441812470,0.28183566809908,
        0.3308819553164567E-01);
      break;
    }
    case 4: {
      base_order=21;
      number_points=1+8*3+10*6; // = 85
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.2761042699769952E-01);
      set1(qp,qw,0.49924967533779,0.1779029547326740E-02);
      set1(qp,qw,0.45293012403052,0.2011239811396117E-01);
      set1(qp,qw,0.39776393795524,0.2681784725933157E-01);
      set1(qp,qw,0.26450020253279,0.2452313380150201E-01);
      set1(qp,qw,0.21101896409208,0.1639457841069539E-01);
      set1(qp,qw,0.10773560717127,0.1479590739864960E-01);
      set1(qp,qw,0.03906908783780,0.4579282277704251E-02);
      set1(qp,qw,0.01117437972933,0.1651826515576217E-02);
      set2(qp,qw,0.00534961818734,0.06354966590835,
        0.2349170908575584E-02);
      set2(qp,qw,0.00795481706620,0.15710691894071,
        0.4465925754181793E-02);
      set2(qp,qw,0.01042239828126,0.39564211436437,
        0.6099566807907972E-02);
      set2(qp,qw,0.01096441479612,0.27316757071291,
        0.6891081327188203E-02);
      set2(qp,qw,0.03856671208546,0.10178538248502,
        0.7997475072478163E-02);
      set2(qp,qw,0.03558050781722,0.44665854917641,
        0.7386134285336024E-02);
      set2(qp,qw,0.04967081636276,0.19901079414950,
        0.1279933187864826E-01);
      set2(qp,qw,0.05851972508433,0.32426118369228,
        0.1725807117569655E-01);
      set2(qp,qw,0.12149778700439,0.20853136321013,
        0.1867294590293547E-01);
      set2(qp,qw,0.14071084494394,0.32317056653626,
        0.2281822405839526E-01);
      break;
    }
    case 5: {
      base_order=26;
      number_points=10*3+16*6; // = 126
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.48602675846341,0.8005581880020417E-02);
      set1(qp,qw,0.43441069933617,0.1594707683239050E-01);
      set1(qp,qw,0.38988913524396,0.1310914123079553E-01);
      set1(qp,qw,0.29844323401980,0.1958300096563562E-01);
      set1(qp,qw,0.23404417233737,0.1647088544153727E-01);
      set1(qp,qw,0.15146833460902,0.8547279074092100E-02);
      set1(qp,qw,0.11273389354599,0.8161885857226492E-02);
      set1(qp,qw,0.07771569209153,0.6121146539983779E-02);
      set1(qp,qw,0.03489309361430,0.2908498264936665E-02);
      set1(qp,qw,0.00725818462093,0.6922752456619963E-03);
      set2(qp,qw,0.00129235270444,0.22721445215336,
        0.1248289199277397E-02);
      set2(qp,qw,0.00539970127212,0.43501055485357,
        0.3404752908803022E-02);
      set2(qp,qw,0.00638400303398,0.32030959927220,
        0.3359654326064051E-02);
      set2(qp,qw,0.00502821150199,0.09175032228001,
        0.1716156539496754E-02);
      set2(qp,qw,0.00682675862178,0.03801083585872,
        0.1480856316715606E-02);
      set2(qp,qw,0.01001619963993,0.15742521848531,
        0.3511312610728685E-02);
      set2(qp,qw,0.02575781317339,0.23988965977853,
        0.7393550149706484E-02);
      set2(qp,qw,0.03022789811992,0.36194311812606,
        0.7983087477376558E-02);
      set2(qp,qw,0.03050499010716,0.08355196095483,
        0.4355962613158041E-02);
      set2(qp,qw,0.04595654736257,0.14844322073242,
        0.7365056701417832E-02);
      set2(qp,qw,0.06744280054028,0.28373970872753,
        0.1096357284641955E-01);
      set2(qp,qw,0.07004509141591,0.40689937511879,
        0.1174996174354112E-01);
      set2(qp,qw,0.08391152464012,0.19411398702489,
        0.1001560071379857E-01);
      set2(qp,qw,0.12037553567715,0.32413434700070,
        0.1330964078762868E-01);
      set2(qp,qw,0.14806689915737,0.22927748355598,
        0.1415444650522614E-01);
      set2(qp,qw,0.19177186586733,0.32561812259598,
        0.1488137956116801E-01);
      break;
    }
    case 6: {
      base_order=31;
      number_points=1+12*3+23*6; // = 175
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.1557996020289920E-01);
      set1(qp,qw,0.49633494178362,0.3177233700534134E-02);
      set1(qp,qw,0.45850216209852,0.1048342663573077E-01);
      set1(qp,qw,0.42450952193729,0.1320945957774363E-01);
      set1(qp,qw,0.38204707005392,0.1497500696627150E-01);
      set1(qp,qw,0.28098784579608,0.1498790444338419E-01);
      set1(qp,qw,0.22734897585403,0.1333886474102166E-01);
      set1(qp,qw,0.17455911150873,0.1088917111390201E-01);
      set1(qp,qw,0.12325842720144,0.8189440660893461E-02);
      set1(qp,qw,0.08008422889220,0.5575387588607785E-02);
      set1(qp,qw,0.04777446740790,0.3191216473411976E-02);
      set1(qp,qw,0.02172051468014,0.1296715144327045E-02);
      set1(qp,qw,0.00476467761544,0.2982628261349172E-03);
      set2(qp,qw,0.00092537119335,0.41529527091331,
        0.9989056850788964E-03);
      set2(qp,qw,0.00138592585556,0.06118990978535,
        0.4628508491732533E-03);
      set2(qp,qw,0.00368241545591,0.16490869013691,
        0.1234451336382413E-02);
      set2(qp,qw,0.00390322342416,0.02503506223200,
        0.5707198522432062E-03);
      set2(qp,qw,0.00323324815501,0.30606446515110,
        0.1126946125877624E-02);
      set2(qp,qw,0.00646743211224,0.10707328373022,
        0.1747866949407337E-02);
      set2(qp,qw,0.00324747549133,0.22995754934558,
        0.1182818815031657E-02);
      set2(qp,qw,0.00867509080675,0.33703663330578,
        0.1990839294675034E-02);
      set2(qp,qw,0.01559702646731,0.05625657618206,
        0.1900412795035980E-02);
      set2(qp,qw,0.01797672125369,0.40245137521240,
        0.4498365808817451E-02);
      set2(qp,qw,0.01712424535389,0.24365470201083,
        0.3478719460274719E-02);
      set2(qp,qw,0.02288340534658,0.16538958561453,
        0.4102399036723953E-02);
      set2(qp,qw,0.03273759728777,0.09930187449585,
        0.4021761549744162E-02);
      set2(qp,qw,0.03382101234234,0.30847833306905,
        0.6033164660795066E-02);
      set2(qp,qw,0.03554761446002,0.46066831859211,
        0.3946290302129598E-02);
      set2(qp,qw,0.05053979030687,0.21881529945393,
        0.6644044537680268E-02);
      set2(qp,qw,0.05701471491573,0.37920955156027,
        0.8254305856078458E-02);
      set2(qp,qw,0.06415280642120,0.14296081941819,
        0.6496056633406411E-02);
      set2(qp,qw,0.08050114828763,0.28373128210592,
        0.9252778144146602E-02);
      set2(qp,qw,0.10436706813453,0.19673744100444,
        0.9164920726294280E-02);
      set2(qp,qw,0.11384489442875,0.35588914121166,
        0.1156952462809767E-01);
      set2(qp,qw,0.14536348771552,0.25981868535191,
        0.1176111646760917E-01);
      set2(qp,qw,0.18994565282198,0.32192318123130,
        0.1382470218216540E-01);
      break;
    }
  }
}
void TriangleWandzuraQuadrature::printOn(ostream &os) const {
  os << "TriangleWandzuraQuadrature: base_order = " << base_order
     << endl;
  Quadrature<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
// David Dunavant,
//  "High Degree Efficient Symmetrical Gaussian Quadrature Rules
//  for the Triangle," 
//  International Journal for Numerical Methods in Engineering,
//  Volume 21, 1985, pages 1129-1148.
TriangleDunavantQuadrature::TriangleDunavantQuadrature(int n) {
  switch (n) {
    case 1: {
      number_points=1;
      base_order=2;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set0(quadrature_points,quadrature_weights,1.);
      break;
    }
    case 2: {
      number_points=3;
      base_order=3;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,1./6.,1./3.);
      break;
    }
    case 3:
    case 4: {
      number_points=2*3; // = 6
      base_order=5;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.091576213509771,0.109951743655322);
      set1(qp,qw,0.445948490915965,0.223381589678011);
      break;
    }
    case 5: {
      number_points=1+2*3; // = 7
      base_order=6;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.225);
      set1(qp,qw,(6.+sqrt(15.))/21.,(155.+sqrt(15.))/1200.);
      set1(qp,qw,(6.-sqrt(15.))/21.,(155.-sqrt(15.))/1200.);
      break;
    }
    case 6: {
      number_points=2*3+6; // = 12
      base_order=7;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.249286745170910,0.116786275726379);
      set1(qp,qw,0.063089014491502,0.050844906370207);
      set2(qp,qw,0.053145049844817,0.310352451033784,
        0.082851075618374);
      break;
    }
    case 7: 
    case 8: {
      number_points=1+3*3+6; // = 16
      base_order=9;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.144315607677787);
      set1(qp,qw,0.459292588292723,0.095091634267285);
      set1(qp,qw,0.170569307751760,0.103217370534718);
      set1(qp,qw,0.050547228317031,0.032458497623198);
      set2(qp,qw,0.008394777409958,0.263112829634638,
        0.027230314174435);
      break;
    }
    case 9: {
      number_points=1+4*3+6; // = 19
      base_order=10;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.097135796282799);
      set1(qp,qw,0.489682519198738,0.031334700227139);
      set1(qp,qw,0.437089591492937,0.077827541004774);
      set1(qp,qw,0.188203535619033,0.079647738927210);
      set1(qp,qw,0.044729513394453,0.025577675658698);
      set2(qp,qw,0.036838412054736,0.221962989160766,
        0.043283539377289);
      break;
    }
    case 10: {
      number_points=1+2*3+3*6; // = 25
      base_order=11;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.090817990382754);
      set1(qp,qw,0.485577633383657,0.036725957756467);
      set1(qp,qw,0.109481575485037,0.045321059435528);
      set2(qp,qw,0.141707219414880,0.307939838764121,
        0.072757916845420);
      set2(qp,qw,0.025003534762686,0.246672560639903,
        0.028327242531057);
      set2(qp,qw,0.009540815400299,0.066803251012200,
        0.009421666963733);
      break;
    }
    case 11:
    case 12: {
      number_points=5*3+3*6; // = 33
      base_order=13;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.488217389773805,0.025731066440455);
      set1(qp,qw,0.439724392294460,0.043692544538038);
      set1(qp,qw,0.271210385012116,0.062858224217885);
      set1(qp,qw,0.127576145541586,0.034796112930709);
      set1(qp,qw,0.021317350453210,0.006166261051559);
      set2(qp,qw,0.115343494534698,0.275713269685514,
        0.040371557766381);
      set2(qp,qw,0.022838332222257,0.281325580989940,
        0.022356773202303);
      set2(qp,qw,0.025734050548330,0.116251915907597,
        0.017316231108659);
      break;
    }
    case 13: {
      number_points=1+6*3+3*6; // = 37
      base_order=14;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.052520923400802);
      set1(qp,qw,0.495048184939705,0.011280145209330);
      set1(qp,qw,0.468716635109574,0.031423518362454);
      set1(qp,qw,0.414521336801277,0.047072502504194);
      set1(qp,qw,0.229399572042831,0.047363586536355);
      set1(qp,qw,0.114424495196330,0.031167529045794);
      set1(qp,qw,0.024811391363459,0.007975771465074);
      set2(qp,qw,0.094853828379579,0.268794997058761,
        0.036848402728732);
      set2(qp,qw,0.018100773278807,0.291730066734288,
        0.017401463303822);
      set2(qp,qw,0.022233076674090,0.126357385491669,
        0.015521786839045);
      break;
    }
    case 14: {
      number_points=6*3+4*6; // = 42
      base_order=15;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.488963910362179,0.021883581369429);
      set1(qp,qw,0.417644719340454,0.032788353544125);
      set1(qp,qw,0.273477528308839,0.051774104507292);
      set1(qp,qw,0.177205532412543,0.042162588736993);
      set1(qp,qw,0.061799883090873,0.014433699669777);
      set1(qp,qw,0.019390961248701,0.004923403602400);
      set2(qp,qw,0.057124757403648,0.172266687821356,
        0.024665753212564);
      set2(qp,qw,0.092916249356972,0.336861459796345,
        0.038571510787061);
      set2(qp,qw,0.014646950055654,0.298372882136258,
        0.014436308113534);
      set2(qp,qw,0.001268330932872,0.118974497696957,
        0.005010228838501);
      break;
    }
    case 15: 
    case 16: 
    case 17: {
      number_points=1+8*3+6*6; // = 61
      base_order=18;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.033437199290803);
      set1(qp,qw,0.497170540556774,0.005093415440507);
      set1(qp,qw,0.482176322624625,0.014670864527638);
      set1(qp,qw,0.450239969020782,0.024350878353672);
      set1(qp,qw,0.400266239377397,0.031107550868969);
      set1(qp,qw,0.252141267970953,0.031257111218620);
      set1(qp,qw,0.162047004658461,0.024815654339665);
      set1(qp,qw,0.075875882260746,0.014056073070557);
      set1(qp,qw,0.015654726967822,0.003194676173779);
      set2(qp,qw,0.010186928826919,0.334319867363658,
        0.008119655318993);
      set2(qp,qw,0.135440871671036,0.292221537796944,
        0.026805742283163);
      set2(qp,qw,0.054423924290583,0.319574885423190,
        0.018459993210822);
      set2(qp,qw,0.012868560833637,0.190704224192292,
        0.008476868534328);
      set2(qp,qw,0.067165782413524,0.180483211648746,
        0.018292796770025);
      set2(qp,qw,0.014663182224828,0.080711313679564,
        0.006665632004165);
      break;
    }
    case 18: 
    case 19: 
    default:
    {
      number_points=1+8*3+8*6; // = 73
      base_order=20;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.032906331388919);
      set1(qp,qw,0.489609987073006,0.010330731891272);
      set1(qp,qw,0.454536892697893,0.022387247263016);
      set1(qp,qw,0.401416680649431,0.030266125869468);
      set1(qp,qw,0.255551654403098,0.030490967802198);
      set1(qp,qw,0.177077942152130,0.024159212741641);
      set1(qp,qw,0.110061053227952,0.016050803586801);
      set1(qp,qw,0.055528624251840,0.008084580261784);
      set1(qp,qw,0.012621863777229,0.002079362027485);
      set2(qp,qw,0.003611417848412,0.395754787356943,
        0.003884876904981);
      set2(qp,qw,0.134466754530780,0.307929983880436,
        0.025574160612022);
      set2(qp,qw,0.014446025776115,0.264566948406520,
        0.008880903573338);
      set2(qp,qw,0.046933578838178,0.358539352205951,
        0.016124546761731);
      set2(qp,qw,0.002861120350567,0.157807405968595,
        0.002491941817491);
      set2(qp,qw,0.223861424097916,0.075050596975911,
        0.018242840118951);
      set2(qp,qw,0.034647074816760,0.142421601113383,
        0.010258563736199);
      set2(qp,qw,0.010161119296278,0.065494628082938,
        0.003799928855302);
      break;
    }
  }
}
void TriangleDunavantQuadrature::printOn(ostream &os) const {
  os << "TriangleDunavantQuadrature: base_order = " << base_order
     << endl;
  Quadrature<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
// Mark Taylor, Beth Wingate, Rachel Vincent,
// An Algorithm for Computing Fekete Points in the Triangle,
// SIAM Journal on Numerical Analysis,
// Volume 38, Number 5, 2000, pages 1707-1720.
TriangleFeketeQuadrature::TriangleFeketeQuadrature(int n) {
  switch (n) {
    case 1: {
      number_points=1+3+6; // = 10
      base_order=4;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.45);
      set1(qp,qw,0.,1./60.);
      set2(qp,qw,0.,(5.-sqrt(5.))/10.,1./12.);
      break;
    }
    case 2: {
      number_points=1+3*3+3*6; // = 28
      base_order=7;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.2178563571);
      set1(qp,qw,0.1063354684,0.5*0.1104193374);
      set1(qp,qw,0.5000000000,0.5*0.0358939762);
      set1(qp,qw,0.0000000000,0.5*0.0004021278);
      set2(qp,qw,0.1171809171,0.3162697959,0.5*0.1771348660);
      set2(qp,qw,0.,0.2655651402,0.5*0.0272344079);
      set2(qp,qw,0.,0.0848854223,0.5*0.0192969460);
      break;
    }
    case 3: {
      number_points=1+4*3+7*6; // = 55
      base_order=10;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.1096011288);
      set1(qp,qw,0.1704318201,0.5*0.0767491008);
      set1(qp,qw,0.4699587644,0.5*0.0646677819);
      set1(qp,qw,0.0489345696,0.5*0.0276211659);
      set1(qp,qw,0.,0.5*0.0013925011);
      set2(qp,qw,0.1784337588,0.3252434900,0.5*0.0933486453);
      set2(qp,qw,0.0588564879,0.3010242110,0.5*0.0619010169);
      set2(qp,qw,0.0551758079,0.1543901944,0.5*0.0437466450);
      set2(qp,qw,0.,0.4173602935,0.5*0.0114553907);
      set2(qp,qw,0.,0.2610371960,0.5*0.0093115568);
      set2(qp,qw,0.,0.1306129092,0.5*0.0078421987);
      set2(qp,qw,0.,0.0402330070,0.5*0.0022457501);
      break;
    }
    case 4: {
      number_points=1+6*3+12*6; // = 91
      base_order=13;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.0626245179);
      set1(qp,qw,0.4005558262,0.5*0.0571359417);
      set1(qp,qw,0.2618405201,0.5*0.0545982307);
      set1(qp,qw,0.0807386775,0.5*0.0172630326);
      set1(qp,qw,0.0336975736,0.5*0.0142519606);
      set1(qp,qw,0.5,0.5*0.0030868485);
      set1(qp,qw,0.,0.5*0.0004270742);
      set2(qp,qw,0.1089969290,0.3837518758,0.5*0.0455876390);
      set2(qp,qw,0.1590834479,0.2454317980,0.5*0.0496701966);
      set2(qp,qw,0.0887037176,0.1697134458,0.5*0.0387998322);
      set2(qp,qw,0.0302317829,0.4071849276,0.5*0.0335323983);
      set2(qp,qw,0.0748751152,0.2874821712,0.5*0.0268431561);
      set2(qp,qw,0.0250122615,0.2489279690,0.5*0.0237377452);
      set2(qp,qw,0.0262645218,0.1206826354,0.5*0.0177255972);
      set2(qp,qw,0.,0.3753565349,0.5*0.0043097313);
      set2(qp,qw,0.,0.2585450895,0.5*0.0028258057);
      set2(qp,qw,0.,0.1569057655,0.5*0.0030994935);
      set2(qp,qw,0.,0.0768262177,0.5*0.0023829062);
      set2(qp,qw,0.,0.0233450767,0.5*0.0009998683);
      break;
    }
    case 5: 
    case 6: {
      number_points=1+7*3+19*6; // = 136
      base_order=16;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.0459710878);
      set1(qp,qw,0.4206960976,0.5*0.0384470625);
      set1(qp,qw,0.2260541354,0.5*0.0386013566);
      set1(qp,qw,0.1186657611,0.5*0.0224308157);
      set1(qp,qw,0.4761452137,0.5*0.0243531004);
      set1(qp,qw,0.0531173538,0.5*0.0094392654);
      set1(qp,qw,0.0219495841,0.5*0.0061105652);
      set1(qp,qw,0.,0.0001283162);
      set2(qp,qw,0.2379370518,0.3270403780,0.5*0.0346650571);
      set2(qp,qw,0.1585345951,0.3013819154,0.5*0.0305412307);
      set2(qp,qw,0.0972525649,0.3853507643,0.5*0.0262101254);
      set2(qp,qw,0.0875150140,0.2749910734,0.5*0.0265367617);
      set2(qp,qw,0.1339547708,0.1975591066,0.5*0.0269859772);
      set2(qp,qw,0.0475622627,0.3524012205,0.5*0.0172635676);
      set2(qp,qw,0.0596194677,0.1978887556,0.5*0.0188795851);
      set2(qp,qw,0.0534939782,0.1162464503,0.5*0.0158224870);
      set2(qp,qw,0.0157189888,0.4176001732,0.5*0.0127170850);
      set2(qp,qw,0.0196887324,0.2844332752,0.5*0.0164489660);
      set2(qp,qw,0.0180698489,0.1759511193,0.5*0.0120018620);
      set2(qp,qw,0.0171941515,0.0816639421,0.5*0.0072268907);
      set2(qp,qw,0.,0.4493368632,0.5*0.0023599161);
      set2(qp,qw,0.,0.3500847655,0.5*0.0017624674);
      set2(qp,qw,0.,0.2569702891,0.5*0.0018648017);
      set2(qp,qw,0.,0.1738056486,0.5*0.0012975716);
      set2(qp,qw,0.,0.1039958541,0.5*0.0018506035);
      set2(qp,qw,0.,0.0503997335,0.5*0.0009919379);
      set2(qp,qw,0.,0.0152159769,0.5*0.0004893506);
      break;
    }
    case 7: 
    {
      number_points=1+9*3+27*6; // = 190
      base_order=19;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<2> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.5*0.0326079297);
      set1(qp,qw,0.4099034502,0.5*0.0288093886);
      set1(qp,qw,0.2438647767,0.5*0.0279490452);
      set1(qp,qw,0.1512564554,0.5*0.0174438045);
      set1(qp,qw,0.4594655253,0.5*0.0203594338);
      set1(qp,qw,0.0832757649,0.5*0.0113349170);
      set1(qp,qw,0.0369065587,0.5*0.0046614185);
      set1(qp,qw,0.0149574850,0.5*0.0030346239);
      set1(qp,qw,0.5,0.5*0.0012508731);
      set1(qp,qw,0.,0.5*0.0000782945);
      set2(qp,qw,0.2515553103,0.3292984162,0.5*0.0255331366);
      set2(qp,qw,0.1821465920,0.3095465041,0.5*0.0235716330);
      set2(qp,qw,0.1246901255,0.3789288931,0.5*0.0206304700);
      set2(qp,qw,0.1179441386,0.2868915642,0.5*0.0204028340);
      set2(qp,qw,0.1639418454,0.2204868669,0.5*0.0215105697);
      set2(qp,qw,0.0742549663,0.3532533654,0.5*0.0183482070);
      set2(qp,qw,0.0937816771,0.2191980979,0.5*0.0174161032);
      set2(qp,qw,0.0890951387,0.1446273457,0.5*0.0155972434);
      set2(qp,qw,0.0409065243,0.4360543636,0.5*0.0119269616);
      set2(qp,qw,0.0488675890,0.2795984854,0.5*0.0147074804);
      set2(qp,qw,0.0460342127,0.2034211147,0.5*0.0116182830);
      set2(qp,qw,0.0420687187,0.1359040280,0.5*0.0087639138);
      set2(qp,qw,0.0116377940,0.4336892286,0.5*0.0098563528);
      set2(qp,qw,0.0299062187,0.3585587824,0.5*0.0096342355);
      set2(qp,qw,0.0132313129,0.2968103667,0.5*0.0086477936);
      set2(qp,qw,0.0136098469,0.2050279257,0.5*0.0083868302);
      set2(qp,qw,0.0124869684,0.1232146223,0.5*0.0062576643);
      set2(qp,qw,0.0365197797,0.0805854893,0.5*0.0077839825);
      set2(qp,qw,0.0118637765,0.0554881302,0.5*0.0031415239);
      set2(qp,qw,0.,0.4154069883,0.5*0.0006513246);
      set2(qp,qw,0.,0.3332475761,0.5*0.0021137942);
      set2(qp,qw,0.,0.2558853572,0.5*0.0004393452);
      set2(qp,qw,0.,0.1855459314,0.5*0.0013662119);
      set2(qp,qw,0.,0.1242528987,0.5*0.0003331251);
      set2(qp,qw,0.,0.0737697111,0.5*0.0011613225);
      set2(qp,qw,0.,0.0355492359,0.5*0.0004342867);
      set2(qp,qw,0.,0.0106941169,0.5*0.0002031499);
      break;
    }
  }
}
void TriangleFeketeQuadrature::printOn(ostream &os) const {
  os << "TriangleFeketeQuadrature: base_order = " << base_order
     << endl;
  Quadrature<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
// Peter Silvester,
// Symmetric Quadrature Formulae for Simplexes,
// Mathematics of Computation,
// Volume 24, Number 109, January 1970, pages 95-100.
TetrahedronNewtonCotesQuadrature::TetrahedronNewtonCotesQuadrature(
int np_1d) {
  switch (np_1d) {
    case 1:
    case 2: { 
      base_order=2;
      number_points=4;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,0.,1./4.);
      break;
    }
    case 4:
    default: {
      base_order=4;
      number_points=2*4;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,0.025);
      set1(qp,qw,1./3.,0.225);
      break;
    }
  }
}
void TetrahedronNewtonCotesQuadrature::printOn(ostream &os) const {
  os << "TetrahedronNewtonCotesQuadrature: base_order = " << base_order
     << endl;
  Quadrature<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
TetrahedronFullNewtonCotesQuadrature::TetrahedronFullNewtonCotesQuadrature(
int np_1d) {
  switch (np_1d) {
    case 1:
    case 2: { 
      base_order=2;
      number_points=4; // vertices
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,0.,1./4.);
      break;
    }
    case 3: { // negative weight
      base_order=3;
      number_points=4+6; // = 10
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,-0.05); // vertices
      setPairs(qp,qw,0.,0.2); // midpoints of edges
      break;
    }
    case 4: {
      base_order=4;
      number_points=2*4+12; // = 20
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,0.025); // vertices
      set1(qp,qw,1./3.,0.225); // lattice points interior to faces
      set2(qp,qw,1./3.,2./3.,0.); // interior lattice points on edges
      break;
    }
    case 5: { // negative weight
      base_order=5;
      number_points=1+4+6+2*12; // = 35
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.3047619047619047894); // interior lattice point
      set1(qp,qw,0.,-0.011904761904761904101); // vertices
      setPairs(qp,qw,0.,-0.028571428571428570536); // midpoints of edges
      set2(qp,qw,0.25,0.75,0.038095238095238098674); // edge lattice points
      set2(qp,qw,0.,0.5,0.038095238095238098674); // interior to faces
      break;
    }
    case 6: { // negative weight
      base_order=6;
      number_points=2*4+4*12; // = 56
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,0.0081845238095238099368);
      set1(qp,qw,0.2,0.093005952380952383596);
      set2(qp,qw,0.2,0.8,-0.008680555555555555941);
      set2(qp,qw,0.4,0.6,0.008680555555555555941);
      set2(qp,qw,0.,0.6,0.068204365079365072977);
      set2(qp,qw,0.,0.2,-0.018601190476190476025);
      break;
    }
    case 7: { // negative weight
      base_order=7;
      number_points=3*4+2*6+3*12+24; // =84
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,-0.005);
      set1(qp,qw,1./6.,0.12857142857142855874);
      set1(qp,qw,1./3.,-0.032142857142857139685);
      setPairs(qp,qw,0.,0.028571428571428570536);
      setPairs(qp,qw,1./6.,0.);
      set2(qp,qw,1./6.,5./6.,0.01714285714285714371);
      set2(qp,qw,1./3.,2./3.,-0.02142857142857142877);
      set2(qp,qw,0.,2./3.,0.);
      set3(qp,qw,0.,1./6.,1./3.,0.02142857142857142877);
      break;
    }
    default:
      OBSOLETE("case not programmed");
  }
}
void TetrahedronFullNewtonCotesQuadrature::printOn(ostream &os) const {
  os << "TetrahedronFullNewtonCotesQuadrature: base_order = " << base_order
     << endl;
  Quadrature<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
TetrahedronGaussianQuadrature::TetrahedronGaussianQuadrature(int ord) { 
  switch (ord) {
    case 1:
    case 2: { // keast 1 = felippa 1
//    exact for | alpha | < 2
      base_order=2;
      number_points=1;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set0(quadrature_points,quadrature_weights,1.);
      break;
    }
    case 3: { // keast 2 = felippa 2
//    exact for | alpha | < 3 and [1,1,1,0]
      base_order=3;
      number_points=4;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,1./(5.+sqrt(5.)),0.25);
      break;
    }
    case 4: { // newton-cotes 4 = felippa 3 
//    exact for | alpha | < 4 
      base_order=4;
      number_points=2*4; // = 8
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,0.025);
      set1(qp,qw,1./3.,0.225);
      break;
    }
    case 5: // keast 5 is unstable, felippa 6 is inefficient
    case 6: { // felippa 5
      base_order=6;
      number_points=2*4+6; // = 14
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.092735250310891226402,0.073493043116361949544);
      set1(qp,qw,0.31088591926330060980,0.11268792571801585080);
      setPairs(qp,qw,0.045503704125649649492,0.042546020777081466438);
      break;
    }
    case 7: { // keast 8 = felippa 9
      base_order=7;
      number_points=3*4+12; // = 24
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.21460287125991520293,0.039922750257869636194);
      set1(qp,qw,0.040673958534611353116,0.010077211055345822612);
      set1(qp,qw,0.32233789014227551034,0.055357181543927398338);
      set2(qp,qw,0.26967233145831580803,0.60300566479164914137,
        0.048214285714285714286);
      break;
    }
  }
}
void TetrahedronGaussianQuadrature::printOn(ostream &os) const {
  os << "TetrahedronGaussianQuadrature: base_order = " << base_order
     << endl;
  Quadrature<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
// mapping is
//   | x_0 | = | xi_0 ( 1 - xi_1 / 2 - xi_2 / 2 + xi_1 xi_2 / 3 ) |
//   | x_1 | = | xi_1 ( 1 - xi_2 / 2 - xi_0 / 2 + xi_2 xi_0 / 3 ) |
//   | x_2 | = | xi_2 ( 1 - xi_0 / 2 - xi_1 / 2 + xi_0 xi_1 / 3 ) |
//mapping determinant is linear, 
//  so Gaussian quad on rectangle is applied to integrand times a linear
//  so order is one less than Gaussian quadrature on a rectangle
TetrahedronProductGaussianQuadrature::TetrahedronProductGaussianQuadrature(
int n) {
  if (n<2) n=2;
  GaussianQuadrature<3> gq(n);
  number_points=gq.numberPoints();
  base_order=gq.order()-2;
  quadrature_points.allocate(number_points);
  quadrature_weights.allocate(number_points);
  for (int p=0;p<number_points;p++) {
    Point<3> &qp=quadrature_points[p];
    const Point<3> &gq_qp=gq.points()[p];
    Tensor<2,3> m;
    m[0][0]=1.-0.5*(gq_qp[1]+gq_qp[2])+gq_qp[1]*gq_qp[2]/3.;
    m[1][0]=gq_qp[1]*(-0.5+gq_qp[2]/3.);
    m[2][0]=gq_qp[2]*(-0.5+gq_qp[1]/3.);
    m[0][1]=gq_qp[0]*(-0.5+gq_qp[2]/3.);
    m[1][1]=1.-0.5*(gq_qp[0]+gq_qp[2])+gq_qp[0]*gq_qp[2]/3.;
    m[2][1]=gq_qp[2]*(-0.5+gq_qp[0]/3.);
    m[0][2]=gq_qp[0]*(-0.5+gq_qp[1]/3.);
    m[1][2]=gq_qp[1]*(-0.5+gq_qp[0]/3.);
    m[2][2]=1.-0.5*(gq_qp[0]+gq_qp[1])+gq_qp[0]*gq_qp[1]/3.;
    REAL det=determinant(m);
    qp[0]=gq_qp[0]*m[0][0];
    qp[1]=gq_qp[1]*m[1][1];
    qp[2]=gq_qp[2]*m[2][2];
    quadrature_weights[p]=gq.weights()[p]*det;
  }
}
void TetrahedronProductGaussianQuadrature::printOn(ostream &os) const {
  os << "TetrahedronProductGaussianQuadrature: base_order = " 
     << base_order << endl;
  Quadrature<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
TetrahedronLobattoQuadrature::TetrahedronLobattoQuadrature(int np_1d) {
  switch (np_1d) {
    case 2: 
    default: {
      base_order=2;
      number_points=4;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      quadrature_points[0]=Point<3>(1.,0.,0.);
      quadrature_points[1]=Point<3>(0.,1.,0.);
      quadrature_points[2]=Point<3>(0.,0.,1.);
      quadrature_points[3]=Point<3>(0.,0.,0.);
      quadrature_weights.initialize(1./24.);
      break;
    }
  }
}
void TetrahedronLobattoQuadrature::printOn(ostream &os) const {
  os << "TetrahedronLobattoQuadrature: base_order = " << base_order
     << endl;
  Quadrature<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
// Patrick Keast,
// Moderate Degree Tetrahedral Quadrature Formulas,
// Computer Methods in Applied Mechanics and Engineering,
// Volume 55, Number 3, May 1986, pages 339-348.
TetrahedronKeastQuadrature::TetrahedronKeastQuadrature(int n) {
  switch (n) {
    case 1: {
      base_order=2;
      number_points=1;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set0(quadrature_points,quadrature_weights,1.);
      break;
    }
    case 2: {
      base_order=3;
      number_points=4;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,(1.-sqrt(.2))*0.25,1./4.);
      break;
    }
    case 3: 
    case 4: {
      base_order=4;
      number_points=4+6; // = 10
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.143856471934385200,6.*0.0362941783134009000);
      setPairs(qp,qw,0.,6.*0.00358165890217718333);
      break;
    }
    case 5: {
      base_order=5;
      number_points=1+4+6; // = 11
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,-6.*0.0131555555555555556);
      set1(qp,qw,0.0714285714285714285,6.*0.00762222222222222222);
      setPairs(qp,qw,0.100596423833200785,6.*0.0248888888888888889);
      break;
    }
    case 6: {
      base_order=5;
      number_points=2*4+6; // = 14
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.100526765225204467,6.*0.0147649707904967828);
      set1(qp,qw,0.314372873493192195,6.*0.0221397911142651221);
      setPairs(qp,qw,0.,6.*0.00317460317460317450);
      break;
    }
    case 7: {
      base_order=6;
      number_points=1+2*4+6; // = 15
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,6.*0.0302836780970891856);
      set1(qp,qw,1./3.,6.*0.00602678571428571597);
      set1(qp,qw,1./11.,6.*0.0116452490860289742);
      setPairs(qp,qw,0.0665501535736642813,6.*0.0109491415613864534);
      break;
    }
    case 8: 
    default: {
      base_order=7;
      number_points=3*4+12; // = 24
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.214602871259151684,6.*0.00665379170969464506);
      set1(qp,qw,0.0406739585346113397,6.*0.00167953517588677620);
      set1(qp,qw,0.322337890142275646,6.*0.00922619692394239843);
      set2(qp,qw,0.269672331458315867,0.603005664791649076,
        6.*0.00803571428571428248);
      break;
    }
  }
}
void TetrahedronKeastQuadrature::printOn(ostream &os) const {
  os << "TetrahedronKeastQuadrature: base_order = " << base_order
     << endl;
  Quadrature<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
// Carlos Felippa,
// A compendium of FEM integration formulas for symbolic work,
// Engineering Computation,
// Volume 21, Number 8, 2004, pages 867-890.
TetrahedronFelippaQuadrature::TetrahedronFelippaQuadrature(int n) {
  switch (n) {
    case 1: {
      base_order=2;
      number_points=1;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set0(quadrature_points,quadrature_weights,1.);
      break;
    }
    case 2: {
      base_order=3;
      number_points=4;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      set1(quadrature_points,quadrature_weights,(1.-sqrt(.2))*0.25,1./4.);
      break;
    }
    case 3: {
      base_order=4;
      number_points=2*4;
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.32805469671142664734,0.13852796651186214232);
      set1(qp,qw,0.10695227393293068277,0.11147203348813785768);
      break;
    }
    case 4: {
      base_order=4;
      number_points=2*4; // = 8
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.,0.025);
      set1(qp,qw,1./3.,0.225);
      break;
    }
    case 5: {
      base_order=6;
      number_points=2*4+6; // = 14
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.092735250310891226402,0.073493043116361949544);
      set1(qp,qw,0.31088591926330060980,0.11268792571801585080);
      setPairs(qp,qw,0.045503704125649649492,0.042546020777081466438);
      break;
    }
    case 6: {
      base_order=5;
      number_points=2*4+6; // = 14
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.31437287349319219275,0.13283874668559071814);
      set1(qp,qw,0.10052676522520447969,0.088589824742980710434);
      setPairs(qp,qw,0.,0.019047619047619047619);
      break;
    }
    case 7: {
      base_order=6;
      number_points=1+2*4+6; // = 15
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.11851851851851851852);
      set1(qp,qw,0.091971078052723032789,0.071937083779018620010);
      set1(qp,qw,0.31979362782962990839,0.069068207226272385281);
      setPairs(qp,qw,0.056350832689629155741,0.052910052910052910053);
      break;
    }
    case 8: {
      base_order=6;
      number_points=1+2*4+6; // = 15
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.18170206858253504234);
      set1(qp,qw,1./3.,0.036160714285714285714);
      set1(qp,qw,1./11.,0.069871494516173816465);
      setPairs(qp,qw,0.066550153573664298240,0.065694849368318756074);
      break;
    }
    case 9: {
      base_order=7;
      number_points=3*4+12; // = 24
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set1(qp,qw,0.21460287125991520293,0.039922750257869636194);
      set1(qp,qw,0.040673958534611353116,0.010077211055345822612);
      set1(qp,qw,0.32233789014227551034,0.055357181543927398338);
      set2(qp,qw,0.26967233145831580803,0.60300566479164914137,
        0.048214285714285714286);
      break;
    }
    default:
      OBSOLETE("case not programmed");
  }
}
void TetrahedronFelippaQuadrature::printOn(ostream &os) const {
  os << "TetrahedronFelippaQuadrature: base_order = " << base_order
     << endl;
  Quadrature<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
// Axel Grundmann, Michael Moeller,
// Invariant Integration Formulas for the N-Simplex 
// by Combinatorial Methods,
// SIAM Journal on Numerical Analysis,
// Volume 15, Number 2, April 1978, pages 282-290.
//all involve negative weights!
TetrahedronGrundmannMoellerQuadrature::
TetrahedronGrundmannMoellerQuadrature(int n) {
  switch (n) {
    case 1: {
      base_order=2;
      number_points=1+4; // = 5
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,-0.8);
      set1(qp,qw,1./6.,0.45);
      break;
    }
    case 2: {
      base_order=2;
      number_points=1+2*4+6; // = 15
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,0.26666666666666666);
      set1(qp,qw,0.125,0.30476190476190473388);
      set1(qp,qw,1./6.,-0.57857142857142851433);
      setPairs(qp,qw,0.125,0.30476190476190473388);
      break;
    }
    case 3: 
    default:
    {
      base_order=3;
      number_points=1+4*4+6+12; // = 35
      quadrature_points.allocate(number_points);
      quadrature_weights.allocate(number_points);
      NumPtr<Point<3> > qp(quadrature_points);
      NumPtr<REAL> qw(quadrature_weights);
      set0(qp,qw,-0.050793650793650793607);
      set1(qp,qw,0.1,0.25834986772486778772);
      set1(qp,qw,0.3,0.25834986772486778772);
      set1(qp,qw,0.125,-0.54179894179894172446);
      set1(qp,qw,1./6.,0.32544642857142858094);
      setPairs(qp,qw,0.125,-0.54179894179894172446);
      set2(qp,qw,0.3,0.5,0.25834986772486778772);
      break;
    }
  }
}
void TetrahedronGrundmannMoellerQuadrature::printOn(ostream &os) const {
  os << "TetrahedronGrundmannMoellerQuadrature: base_order = " << base_order
     << endl;
  Quadrature<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// explicit instantiations; note: we need them all for all dimensions
//template class ProjectedQuadrature<1>;
#include <NumPtr.C>
INSTANTIATE_NUMPTR(Point<0>);
INSTANTIATE_NUMPTR(Point<1>);
#if (SPACEDIM>1)
INSTANTIATE_NUMPTR(Point<2>);
#endif
#if (SPACEDIM==3)
INSTANTIATE_NUMPTR(Point<3>);
#endif
