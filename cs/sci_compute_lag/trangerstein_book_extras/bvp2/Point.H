#ifndef POINT_H
#define POINT_H "$Header:$"
//----------------------------  point.h  ---------------------------
//    $Id: point.h,v 1.30 2003/01/03 22:21:26 deal Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  point.h  ---------------------------
//
//modified from deal.II/base/include/base/point.h
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

#include <Const.H>
#include <Errors.H>
#include <Types.H>
#include <cmath>

//compare to DEAL Point
template <int dim> class Point {
  private:
    double data[dim];
  public:
    Point() { for (int i=0;i<dim;i++) data[i]=UNDEFINED; }
    Point (const Point<dim> &p) { 
      for (int i=0;i<dim;i++) data[i]=p.data[i]; 
    }
    explicit Point(double x) { CHECK_TEST(dim==1); data[0]=x; } 
    explicit Point(double x,double y) { 
      CHECK_TEST(dim==2); data[0]=x; data[1]=y; 
    }
    explicit Point(double x,double y,double z) {
      CHECK_TEST(dim==3); data[0]=x; data[1]=y; data[2]=z;
    }
    explicit Point(double w,double x,double y,double z) {
      CHECK_TEST(dim==4); data[0]=w; data[1]=x; data[2]=y; data[3]=z;
    }
    bool operator==(const Point<dim> &p) const {
      for (int i=0;i<dim;i++) {
        if (data[i]!=p.data[i]) return false;
      }
      return true;
    }
    Point<dim>& operator=(const Point<dim> &p) {
      for (int i=0;i<dim;i++) data[i]=p.data[i]; return *this;
    }
    Point<dim>& operator=(double p) {
      for (int i=0;i<dim;i++) data[i]=p; return *this;
    }
    double operator[](int i) const { 
      CHECK_BOUNDS(i,0,dim);
      return data[i]; 
    }
    double& operator[] (int i) { 
      CHECK_BOUNDS(i,0,dim);
      return data[i]; 
    }
    Point<dim>& operator+=(const Point<dim> &p) {
      for (int i=0;i<dim;i++) data[i]+=p.data[i];
      return *this;
    }
    Point<dim> operator+(const Point<dim> &p) const {
      Point<dim> q(*this); q+=p; return q; 
    }
    Point<dim>& operator-=(const Point<dim> &p) {
      for (int i=0;i<dim;i++) data[i]-=p.data[i];
      return *this;
    }
    Point<dim> operator-(const Point<dim> &p) const {
      Point<dim> q(*this); q-=p; return q; 
    }
    Point<dim> operator-() const {
      Point<dim> q(*this); 
      for (int i=0;i<dim;i++) q.data[i]=-q.data[i];
      return q;
    }
    Point<dim>& operator*=(double x) {
      for (int i=0;i<dim;i++) data[i]*=x;
      return *this;
    }
    Point<dim> operator*(double x) const {
      Point<dim> q(*this); q*=x; return q;
    }
    double operator*(const Point<dim> &p) const {
      double prod=0.; 
      for (int i=0;i<dim;i++) prod+=data[i]*p.data[i];
      return prod;
    }
    Point<dim> operator/=(double x) {
      CHECK_TEST(abs(x)!=0.);
      for (int i=0;i<dim;i++) data[i]/=x;
      return *this;
    }
    Point<dim> operator/(double x) const {
      Point<dim> q(*this); q/=x; return q;
    }
    void initialize(double r) { 
      for (int i=0;i<dim;++i) data[i]=r;
    }
    double square() const { return (*this) * (*this); }
    double distance(const Point<dim> &p) const {
      return sqrt((*this-p).square());
    }
};

template <int dim> inline ostream& operator<<(ostream &os,
const Point<dim> &p) {
  for (int i=0;i<dim-1;++i) os << p[i] << ' ';
  os << p[dim-1];
  return os;
}

template <> class Point<0> {
  private:
    double bogus;
  public:
    Point() : bogus(UNDEFINED) {}
    Point (const Point<0> &p) : bogus(UNDEFINED) { }
    Point(double) { OBSOLETE("incompatible dim"); }
    Point(double,double) { OBSOLETE("incompatible dim"); }
    Point(double,double,double) { OBSOLETE("incompatible dim"); }
    bool operator==(const Point<0> &p) const { return true; }
    double operator[](int) const { return bogus; }
    double& operator[] (int) { return bogus; }
    Point<0>& operator+=(const Point<0>&) { return *this; }
    Point<0> operator+(const Point<0>&) const { return *this; }
    Point<0>& operator-=(const Point<0> &p) { return *this; }
    Point<0> operator-(const Point<0> &p) const { return *this; }
    Point<0> operator-() const { return *this; }
    Point<0>& operator*=(double) { return *this; }
    Point<0> operator*(double) const { return *this; }
    double operator*(const Point<0> &p) const { return bogus; }
    Point<0> operator/=(double) { return *this; }
    Point<0> operator/(double) const { return *this; }
    double square() const { return bogus; }
    double distance(const Point<0> &p) const { return bogus; }
};

template <> inline ostream& operator<<(ostream &os,const Point<0> &p) {
  return os;
}

#endif
