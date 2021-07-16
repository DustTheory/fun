#ifndef QUADRATURE_H
#define QUADRATURE_H "$Header:$"
//----------------------------  quadrature.h  ---------------------------
//    $Id: quadrature.h,v 1.49 2003/05/16 21:11:23 wolf Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  quadrature.h  ---------------------------
//----------------------------  quadrature_lib.h  -----------------------
//    $Id: quadrature_lib.h,v 1.25 2003/02/17 08:16:15 hartmann Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  quadrature_lib.h  -----------------------
//
//modified from deal.II/base/include/base/quadrature.h and quadrature_lib.h
// by John Trangenstein, August 2009
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

#include <NumPtr.H>
#include <Point.H>

//compare to DEAL Quadrature
template<int dim> class Quadrature {
  protected:
    NumPtr<Point<dim> > quadrature_points;
    NumPtr<double> quadrature_weights;
    int number_points;
  public:
    Quadrature() : number_points(0) {}
    Quadrature(int n) : number_points(n),quadrature_points(n),
    quadrature_weights(n) {
      quadrature_weights.initialize(UNDEFINED);
    }

    Quadrature(const Quadrature<dim-1>&, const Quadrature<1>&); 
    Quadrature(const NumPtr<Point<dim> > &p,const NumPtr<double> &w) :
    number_points(p.getNumber()),quadrature_points(p.getNumber()),
    quadrature_weights(w.getNumber()) {
      CHECK_SAME(number_points,w.getNumber());
      quadrature_points.copyFrom(p);
      quadrature_weights.copyFrom(w);
    }
    Quadrature(const NumPtr<Point<dim> > &p) : 
    number_points(p.getNumber()),quadrature_points(p.getNumber()),
    quadrature_weights(p.getNumber()) {
      quadrature_points.copyFrom(p);
      quadrature_weights.initialize(UNDEFINED); 
    }
    Quadrature(const Point<dim> &p) : number_points(1),
    quadrature_points(1),quadrature_weights(1) {
      quadrature_points[0]=p; quadrature_weights[0]=1.;
    }
    virtual ~Quadrature() {;}
    Quadrature& operator=(const Quadrature &q) { 
      number_points=q.number_points;
      quadrature_points.cleanup();
      quadrature_points.allocate(q.quadrature_points.getNumber());
      quadrature_points.copyFrom(q.quadrature_points);
      quadrature_weights.cleanup();
      quadrature_weights.allocate(q.quadrature_weights.getNumber());
      quadrature_weights.copyFrom(q.quadrature_weights);
      return *this;
    }
    int numberPoints() const { return number_points; }
    const Point<dim>& point(int i) const {
      return quadrature_points[i];
    }
    const NumPtr<Point<dim> >& points() const {
      return quadrature_points;
    }
    void setPoint(int i,const Point<dim> &p) { 
      quadrature_points[i]=p; 
    }
    double weight(int i) const { return quadrature_weights[i]; }
    const NumPtr<double>& weights() const { return quadrature_weights; }
//  1+degree of poly integrated exactly:
    virtual int order() const { return 0; };
//  choose constructor argument given desired quadrature order of accuracy
//  quadrature will be ``exact'' for polynomials of degree order-1
    virtual void printOn(ostream&) const;
};
template <> Quadrature<0>::Quadrature();
template <> Quadrature<0>::Quadrature(const Quadrature<-1>&,
                                      const Quadrature<1>&);
template <> Quadrature<1>::Quadrature(const Quadrature<0>&,
                                      const Quadrature<1>&);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int dimensionalRoot(int,int);
#if (SPACEDIM>=2)
int triangleRoot(int);
#endif
#if (SPACEDIM==3)
int tetrahedronRoot(int);
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//compare to DEAL QGauss2, QGauss3, QGauss4, QGauss5, QGauss6, QGauss7
template <int dim> class GaussianQuadrature : public Quadrature<dim> {
  private:
    int np_1d;
  public:
    GaussianQuadrature(int);
    ~GaussianQuadrature() {;}
    int order() const { return 2*np_1d; }
    static int argGivenOrder(int ord) { return (max(1,ord)+1)/2; }
    static int nPointsGivenOrder(int ord) { 
      return __cmath_power(argGivenOrder(ord),dim); 
    }
    static int argGivenNPoints(int n) { return dimensionalRoot(dim,n); }
    virtual void printOn(ostream&) const;
};
//midpoint rule is GaussianQuadrature(1)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//not in DEAL
template <int dim> class LobattoQuadrature : public Quadrature<dim> {
  private:
    int np_1d;
  public:
    LobattoQuadrature(int);
    ~LobattoQuadrature() {;}
    int order() const { return 2*(np_1d-1); }
    static int argGivenOrder(int ord) { return 1+(max(1,ord)+1)/2; }
    static int nPointsGivenOrder(int ord) { 
      return __cmath_power(argGivenOrder(ord),dim); 
    }
    static int argGivenNPoints(int n) { return dimensionalRoot(dim,n); }
    virtual void printOn(ostream&) const;
};
//trapezoidal rule is LobattoQuadrature(1)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//compare to DEAL QTrapez, QMilne, QWeddle
template <int dim> class NewtonCotesQuadrature: public Quadrature<dim> {
  private:
    int np_1d;
  public:
    NewtonCotesQuadrature(int);
    ~NewtonCotesQuadrature() {;}
    int order() const { return ((np_1d+1)/2)*2; }
    static int argGivenOrder(int ord) { 
      return min(8,max(2,((ord+1)/2)*2-1)); 
    }
    static int nPointsGivenOrder(int ord) { 
      return __cmath_power(argGivenOrder(ord),dim); 
    }
    static int argGivenNPoints(int n) { return dimensionalRoot(dim,n); }
    virtual void printOn(ostream&) const;
};
//trapezoidal rule is NewtonCotesQuadrature(2)
//Simpson's rule is NewtonCotesQuadrature(3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//not in DEAL
template <int dim> class ClenshawCurtisQuadrature :
public Quadrature<dim> {
  private:
    int np_1d;
  public:
    ClenshawCurtisQuadrature(int);
    ~ClenshawCurtisQuadrature() {;}
    int order() const { return 2*((np_1d+1)/2); }
    static int argGivenOrder(int ord) { return ord; }
    static int nPointsGivenOrder(int ord) { 
      return __cmath_power(argGivenOrder(ord),dim); 
    }
    static int argGivenNPoints(int n) { return dimensionalRoot(dim,n); }
    virtual void printOn(ostream&) const;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
//compare to DEAL QIterated
class ProductQuadrature : public Quadrature<3> {
  private:
    int base_order;
  public:
    ProductQuadrature(const Quadrature<2>&,const Quadrature<1>&);
    ~ProductQuadrature() {;}
    int order() const { return base_order; }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//see http://www.cs.cornell.edu/home/vavasis/quad.html
//for references on quadrature rules
#if (SPACEDIM>1)
//not in DEAL
class TriangleNewtonCotesQuadrature : public Quadrature<2> {
  private:
    int base_order; 
  public:
    TriangleNewtonCotesQuadrature(int);
    ~TriangleNewtonCotesQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=2) return 2;
      if (ord==3) return 3;
      if (ord>=4) return 5;
    }
    static int nPointsGivenOrder(int ord) { 
      if (ord<=3) return 3;
      if (ord>=4) return 9;
    }
    static int argGivenNPoints(int n) { 
      if (n<=3) return 3;
      if (n<=4) return 4;
      if (n<=9) return 5;
      return 6;
    }
    virtual void printOn(ostream&) const;
};
//like TriangleNewtonCotesQuadrature, but some points have zero weight
//  number of quadrature points = number of equallySpacedLatticePoints
//not in DEAL
class TriangleFullNewtonCotesQuadrature : public Quadrature<2> {
  private:
    int base_order; 
  public:
    TriangleFullNewtonCotesQuadrature(int);
    ~TriangleFullNewtonCotesQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=2) return 2;
      if (ord==3) return 3;
      if (ord==4) return 5;
      if (ord>=5) return 8;
    }
    static int nPointsGivenOrder(int ord) { 
      if (ord<=2) return 3;
      if (ord==3) return 6;
      if (ord==4) return 15;
      return 36;
    }
    static int argGivenNPoints(int n) { return triangleRoot(n); }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
//not in DEAL
class TriangleGaussianQuadrature : public Quadrature<2> {
  private:
    int base_order; 
  public:
    TriangleGaussianQuadrature(int);
    ~TriangleGaussianQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=2) return 2;
      if (ord==3) return 3;
      if (ord<=5) return 5;
      if (ord==6) return 6;
      if (ord<=9) return 9;
      if (ord==10) return 10;
      if (ord==11) return 11;
      if (ord<=13) return 13;
      if (ord==14) return 14;
      if (ord==15) return 15;
      if (ord==16) return 16;
      if (ord<=18) return 18;
      if (ord<=20) return 20;
      if (ord==21) return 21;
      if (ord<=26) return 26;
      return 31;
    }
    static int nPointsGivenOrder(int ord) { 
      if (ord<=2) return 1;
      if (ord==3) return 3;
      if (ord<=5) return 6;
      if (ord==6) return 7;
      if (ord<=9) return 16;
      if (ord==10) return 19;
      if (ord==11) return 25;
      if (ord<=13) return 33;
      if (ord==14) return 37;
      if (ord==15) return 42;
      if (ord==16) return 54;
      if (ord<=18) return 61;
      if (ord<=20) return 73;
      if (ord==21) return 85;
      if (ord<=26) return 126;
      return 175;
    }
    static int argGivenNPoints(int n) { 
      if (n<=1) return 2;
      if (n<=3) return 3;
      if (n<=6) return 5;
      if (n<=7) return 6;
      if (n<=16) return 9;
      if (n<=19) return 10;
      if (n<=25) return 11;
      if (n<=33) return 13;
      if (n<=37) return 14;
      if (n<=42) return 15;
      if (n<=54) return 16;
      if (n<=61) return 18;
      if (n<=73) return 20;
      if (n<=85) return 21;
      if (n<=126) return 26;
      return 31;
    }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
//product Gaussian rule mapped to a triangle
//not in DEAL
class TriangleProductGaussianQuadrature : public Quadrature<2> {
  private:
    int base_order; 
  public:
    TriangleProductGaussianQuadrature(int);
    ~TriangleProductGaussianQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { return max(1,ord/2+1); }
    static int nPointsGivenOrder(int ord) {
      __cmath_power(argGivenOrder(ord),2);
    }
    static int argGivenNPoints(int n) { return dimensionalRoot(2,n); }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
//not in DEAL
class TriangleLobattoQuadrature : public Quadrature<2> {
  private:
    int base_order; 
  public:
    TriangleLobattoQuadrature(int);
    ~TriangleLobattoQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=2) return 2;
      if (ord<=4) return ord;
      if (ord<=7) return 7;
      if (ord<=10) return 10;
      if (ord<=13) return 13;
      if (ord<=16) return 16;
      return 19;
    }
    static int nPointsGivenOrder(int ord) {
      if (ord<=2) return 3;
      if (ord==3) return 4;
      if (ord==4) return 10;
      if (ord<=7) return 28;
      if (ord<=10) return 55;
      if (ord<=13) return 91;
      if (ord<=16) return 136;
      return 19;
    }
    static int argGivenNPoints(int n) { 
      if (n<=3) return 2;
      if (n<=4) return 3;
      if (n<=10) return 4;
      if (n<=28) return 7;
      if (n<=55) return 10;
      if (n<=91) return 13;
      if (n<=136) return 16;
      return 19;
    }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
//not in DEAL
class TriangleWandzuraQuadrature : public Quadrature<2> {
  private:
    int base_order; 
  public:
    TriangleWandzuraQuadrature(int);
    ~TriangleWandzuraQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      return max(1,min(6,(ord+3)/5)); 
    }
    static int nPointsGivenOrder(int ord) {
      if (ord<=6) return 7;
      if (ord<=11) return 25;
      if (ord<=16) return 54;
      if (ord<=21) return 85;
      if (ord<=26) return 126;
      if (ord<=31) return 175;
    }
    static int argGivenNPoints(int n) { 
      if (n<=3) return 1;
      if (n<=25) return 2;
      if (n<=54) return 3;
      if (n<=85) return 4;
      if (n<=126) return 5;
      return 6;
    }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
//not in DEAL
class TriangleDunavantQuadrature : public Quadrature<2> {
  private:
    int base_order; 
  public:
    TriangleDunavantQuadrature(int);
    ~TriangleDunavantQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=2) return 1;
      if (ord==3) return 2;
      if (ord<=5) return 4;
      if (ord==6) return 5;
      if (ord==7) return 6;
      if (ord<=9) return 8;
      if (ord==10) return 9;
      if (ord==11) return 10;
      if (ord<=13) return 12;
      if (ord==14) return 13;
      if (ord==15) return 14;
      if (ord<=18) return 17;
      return 19;
    }
    static int nPointsGivenOrder(int ord) {
      if (ord<=2) return 1;
      if (ord==3) return 3;
      if (ord<=5) return 6;
      if (ord==6) return 7;
      if (ord==7) return 12;
      if (ord<=9) return 16;
      if (ord==10) return 19;
      if (ord==11) return 25;
      if (ord<=13) return 33;
      if (ord==14) return 37;
      if (ord==15) return 42;
      if (ord<=18) return 61;
      return 73;
    }
    static int argGivenNPoints(int n) { 
      if (n<=1) return 1;
      if (n<=3) return 2;
      if (n<=6) return 4;
      if (n<=7) return 5;
      if (n<=12) return 6;
      if (n<=16) return 8;
      if (n<=19) return 9;
      if (n<=25) return 10;
      if (n<=33) return 12;
      if (n<=37) return 13;
      if (n<=42) return 14;
      if (n<=61) return 17;
      return 19;
    }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
//not in DEAL
class TriangleFeketeQuadrature : public Quadrature<2> {
  private:
    int base_order; 
  public:
    TriangleFeketeQuadrature(int);
    ~TriangleFeketeQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=13) return max(1,(ord-1)/3); 
      if (ord<=16) return 6; // rule 5 is unstable
      return 7;
    }
    static int nPointsGivenOrder(int ord) {
      if (ord<=4) return 10;
      if (ord<=7) return 28;
      if (ord<=10) return 55;
      if (ord<=13) return 91;
      if (ord<=16) return 136;
      if (ord<=16) return 190; // rule 5 is unstable
      return 190;
    }
    static int argGivenNPoints(int n) { 
      if (n<=10) return 1;
      if (n<=28) return 2;
      if (n<=55) return 3;
      if (n<=91) return 4;
      if (n<=136) return 6;
      return 7;
    }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
//not in DEAL
class TetrahedronNewtonCotesQuadrature : public Quadrature<3> {
  private:
    int base_order; 
  public:
    TetrahedronNewtonCotesQuadrature(int);
    ~TetrahedronNewtonCotesQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { return (ord<=2 ? 2 : 4); }
    static int nPointsGivenOrder(int ord) { 
      if (ord<=2) return 4;
      return 8;
    }
    static int argGivenNPoints(int n) { 
      if (n<=4) return 2;
      return 4;
    }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
//not in DEAL
class TetrahedronFullNewtonCotesQuadrature : public Quadrature<3> {
  private:
    int base_order; 
  public:
    TetrahedronFullNewtonCotesQuadrature(int);
    ~TetrahedronFullNewtonCotesQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { return max(2,min(7,ord)); } 
    static int nPointsGivenOrder(int ord) { 
      if (ord<=2) return 4;
      if (ord==3) return 10;
      if (ord==4) return 20;
      if (ord==5) return 35;
      if (ord==6) return 56;
      return 84;
    }
    static int argGivenNPoints(int n) { return tetrahedronRoot(n); } 
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
//not in DEAL
class TetrahedronGaussianQuadrature : public Quadrature<3> {
  private:
    int base_order; 
  public:
    TetrahedronGaussianQuadrature(int);
    ~TetrahedronGaussianQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=2) return 2;
      if (ord<=4) return ord;
      if (ord<=6) return 6;
      return 7;
    }
    static int nPointsGivenOrder(int ord) { 
      if (ord<=2) return 1;
      if (ord==3) return 4;
      if (ord==4) return 8;
      if (ord<=6) return 14;
      return 24;
    }
    static int argGivenNPoints(int n) {
      if (n<=1) return 2;
      if (n<=4) return 3;
      if (n<=8) return 4;
      if (n<=14) return 6;
      return 7;
    }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
//product Gaussian rule mapped to a tetrahedron
//not in DEAL
class TetrahedronProductGaussianQuadrature : public Quadrature<3> {
  private:
    int base_order; 
  public:
    TetrahedronProductGaussianQuadrature(int);
    ~TetrahedronProductGaussianQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { return max(1,(ord+1)/2+1); }
    static int nPointsGivenOrder(int ord) {
      __cmath_power(argGivenOrder(ord),3);
    }
    static int argGivenNPoints(int n) { return dimensionalRoot(3,n); }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
//not in DEAL
class TetrahedronLobattoQuadrature : public Quadrature<3> {
  private:
    int base_order; 
  public:
    TetrahedronLobattoQuadrature(int);
    ~TetrahedronLobattoQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { return 2; }
    static int nPointsGivenOrder(int ord) { return 4; }
    static int argGivenNPoints(int n) { return 2; }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
//not in DEAL
class TetrahedronKeastQuadrature : public Quadrature<3> {
  private:
    int base_order; 
  public:
    TetrahedronKeastQuadrature(int);
    ~TetrahedronKeastQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=2) return 1;
      if (ord==3) return 2;
      if (ord<=5) return ord;
      return min(8,ord+1);
    }
    static int nPointsGivenOrder(int ord) { 
      if (ord<=2) return 1;
      if (ord==3) return 4;
      if (ord==4) return 10;
      if (ord==5) return 11;
      if (ord==6) return 15;
      return 24;
    }
    static int argGivenNPoints(int n) {
      if (n<=1) return 1;
      if (n<=4) return 2;
      if (n<=10) return 4;
      if (n<=11) return 5;
      if (n<=14) return 6;
      if (n<=15) return 7;
      return 8;
    }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
//not in DEAL
class TetrahedronFelippaQuadrature : public Quadrature<3> {
  private:
    int base_order; 
  public:
    TetrahedronFelippaQuadrature(int);
    ~TetrahedronFelippaQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=2) return 1;
      if (ord==3) return 2;
      if (ord==4) return 3;
      if (ord<=6) return 5;
      return 9;
    }
    static int nPointsGivenOrder(int ord) { 
      if (ord<=2) return 1;
      if (ord==3) return 4;
      if (ord==4) return 8;
      if (ord<=6) return 14;
      return 24;
    }
    static int argGivenNPoints(int n) {
      if (n<=1) return 1;
      if (n<=4) return 2;
      if (n<=8) return 3;
      if (n<=14) return 5;
      if (n<=15) return 7;
      return 9;
    }
    virtual void printOn(ostream&) const;
};
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
//not in DEAL
class TetrahedronGrundmannMoellerQuadrature : public Quadrature<3> {
  private:
    int base_order; 
  public:
    TetrahedronGrundmannMoellerQuadrature(int);
    ~TetrahedronGrundmannMoellerQuadrature() {;}
    int order() const { return base_order; }
    static int argGivenOrder(int ord) { 
      if (ord<=2) return 1;
      return 3;
    }
    static int nPointsGivenOrder(int ord) { 
      if (ord<=2) return 5;
      return 35;
    }
    static int argGivenNPoints(int n) {
      if (n<=5) return 1;
      if (n<=15) return 2;
      return 3;
    }
    virtual void printOn(ostream&) const;
};
#endif

#endif
