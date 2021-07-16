#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H "$Header:$"

//----------------------------  polynomial.h  ---------------------------
//    $Id: polynomial.h,v 1.33 2003/05/12 21:45:05 wolf Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  polynomial.h  ---------------------------
//
//modified from deal.II/base/include/base/polynomial.h
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

#include <NumPtr.H>
#include <Types.H>

//modified substantially from DEAL namespace Polynomials
class Polynomial {
  public:
    Polynomial() {;}
    virtual ~Polynomial() {;}
    virtual void values(double,NumPtr<double>&) const=0;
    virtual void slopes(double,NumPtr<double>&) const=0;
    virtual void slope2s(double,NumPtr<double>&) const=0;
    virtual void derivatives(double,int,NumPtr<double>&) const=0;
    virtual void printOn(ostream&) const;
};

//not in DEAL:
class Monomial : public Polynomial {
  private:
    int order;
  public:
    Monomial(int ord) : order(ord) {;}
    virtual ~Monomial() {;}
    int getOrder() const { return order; }
    virtual void values(double,NumPtr<double>&) const;
    virtual void slopes(double,NumPtr<double>&) const;
    virtual void slope2s(double,NumPtr<double>&) const;
    virtual void derivatives(double,int,NumPtr<double>&) const;
    virtual void printOn(ostream&) const;
};

//arbitrary order
//compare to DEAL LagrangeEquidistant, which is limited to degree <= 10
class C0LagrangePolynomial: public Polynomial {
  private:
    int order;
    NumPtr<double> zeros;
    NumPtr<double> den;
  public:
    C0LagrangePolynomial(int);
    C0LagrangePolynomial(const NumPtr<double>&);
    virtual ~C0LagrangePolynomial() {;}
    const NumPtr<double>& zerosValues() const { return zeros; }
    inline static int maxOrder() { return INT_MAX; }
    int getOrder() const { return order; }
    void values(double,NumPtr<double>&) const;
    void slopes(double,NumPtr<double>&) const;
    void slope2s(double,NumPtr<double>&) const;
    void derivatives(double,int,NumPtr<double>&) const;
    void printOn(ostream&) const;
};

//not in DEAL:
class C1LagrangePolynomial: public Polynomial {
  private:
    int order;
    int n;
    NumPtr<double> zeros;
    NumPtr<double> a;
    NumPtr<double> b;
    NumPtr<double> den;
  public:
    C1LagrangePolynomial(int);
    virtual ~C1LagrangePolynomial() {;}
    const NumPtr<double>& zerosValues() const { return zeros; }
    inline static int maxOrder() { return INT_MAX; }
    int getOrder() const { return order; }
    void values(double,NumPtr<double>&) const;
    void slopes(double,NumPtr<double>&) const;
    void slope2s(double,NumPtr<double>&) const;
    void derivatives(double,int,NumPtr<double>&) const;
    void printOn(ostream&) const;
};

//not in DEAL:
class C2LagrangePolynomial: public Polynomial {
  private:
    int order;
    int n;
    NumPtr<double> zeros;
    NumPtr<double> a0;
    NumPtr<double> b0;
    NumPtr<double> c0;
    NumPtr<double> a1;
    NumPtr<double> b1;
    NumPtr<double> den;
  public:
    C2LagrangePolynomial(int);
    virtual ~C2LagrangePolynomial() {;}
    const NumPtr<double>& zerosValues() const { return zeros; }
    inline static int maxOrder() { return INT_MAX; }
    int getOrder() const { return order; }
    void values(double,NumPtr<double>&) const;
    void slopes(double,NumPtr<double>&) const;
    void slope2s(double,NumPtr<double>&) const;
    void derivatives(double,int,NumPtr<double>&) const;
    void printOn(ostream&) const;
};

//compare to DEAL Legendre
class LegendrePolynomial : public Polynomial {
  public:
    LegendrePolynomial() {;}
    virtual ~LegendrePolynomial() {;}
    inline static int maxOrder() { return INT_MAX; }
    void values(double,NumPtr<double>&) const;
    void slopes(double,NumPtr<double>&) const;
    void slope2s(double,NumPtr<double>&) const;
    void derivatives(double,int,NumPtr<double>&) const;
    void printOn(ostream&) const;
};

//compare to DEAL Hierarchical
class HierarchicalPolynomial : public Polynomial {
  public:
    HierarchicalPolynomial() {;}
    virtual ~HierarchicalPolynomial() {;}
    inline static int maxOrder() { return INT_MAX; }
    void values(double,NumPtr<double>&) const;
    void slopes(double,NumPtr<double>&) const;
    void slope2s(double,NumPtr<double>&) const;
    void derivatives(double,int,NumPtr<double>&) const;
    void printOn(ostream&) const;
};  

//not in DEAL:
class HierarchicalSidePolynomial : public Polynomial {
  public:
    HierarchicalSidePolynomial() {;}
    virtual ~HierarchicalSidePolynomial() {;}
    inline static int maxOrder() { return INT_MAX; }
    void values(double,NumPtr<double>&) const;
    void slopes(double,NumPtr<double>&) const;
    void slope2s(double,NumPtr<double>&) const;
    void derivatives(double,int,NumPtr<double>&) const;
    void printOn(ostream&) const;
};  

#endif
