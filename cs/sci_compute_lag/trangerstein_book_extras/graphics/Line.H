//**********************************************************************
// Copyright 2006 John A. Trangenstein
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
#ifndef LINE_H
#define LINE_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/Line.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#if (SPACEDIM>1)
#include <iostream>
#include "Vector3.H"

class Line {
  private:
    Vector3 position;
    Vector3 direction;
  public:
    inline Line(const Vector3 &p,const Vector3 &d) : position(p),
    direction(d) {
      double length=norm(direction);
      CHECK_POSITIVE(length)
      length=ONE/length;
      direction *= length;
    }
    inline virtual ~Line() {;}
    inline Vector3 getPosition() const { return position; }
    inline Vector3 getDirection() const { return direction; }
    double nearestPositionTo(const Line &line) {
      Vector3 proj=
        direction-line.direction*dot(direction,line.direction);
      Vector3 residual=line.position-position;
      double xi1=dot(line.direction,residual);
      residual -= line.direction*xi1;

      double proj_norm=sumSquares(proj);
      if (proj_norm>ZERO) {
        double xi2=dot(proj,residual)/proj_norm;
        position+=direction*xi2;
        residual-=proj*dot(proj,residual)/proj_norm;
      }
      return norm(residual);
    }
    virtual void printOn(ostream &os) const {
      os << "\tLine:" << endl;
    }
};
#endif
#endif
