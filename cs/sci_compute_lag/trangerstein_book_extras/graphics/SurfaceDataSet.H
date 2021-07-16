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
#ifndef SURFACE_DATASET_H
#define SURFACE_DATASET_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/SurfaceDataSet.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"

#if (SPACEDIM==3)
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>
#include "Vector3.H"
#include "VolGraphTool.H"

class VolGraphTool;
class SurfaceDataSet {
  private:
    VolGraphTool *graph_tool;
    const Vector3 &low;
    const Vector3 &high;
    Vector3 center;
    double radius,var_min,var_max;
    int nptsf,nptsc;
    double (*surface_func)(const Vector3&);

  public:
    SurfaceDataSet(VolGraphTool*,const Vector3&,
      const Vector3&,int,int,double (*)(const Vector3&),double,double);
    virtual ~SurfaceDataSet() {}

    double getVarMin() const { return var_min; }
    double getVarMax() const { return var_max; }
    VolGraphTool* getGraphTool() const { return graph_tool; }
    bool drawSurface() const { return true; }

    void plot(bool,const AxisClipPlane*);
    void plot(bool);
    void plot(bool,double);
    virtual void printOn(ostream&) const;
};

#endif
#endif
