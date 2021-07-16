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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/XYGraphTool.C,v 1.1 2009/08/20 17:31:47 johnt Exp $"
#include "XYGraphTool.H"
#include <algorithm>
#include <string.h>
#include <stdlib.h>
#include "Const.H"
//#include "Tracer.H"
# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool::XYGraphTool(const char *str,const char *xlabel,
const char *ylabel,double xxlo,double xxhi,double yylo,double yyhi,
VirtualColormap *cm,char *dn,double maxwinsize,bool mono) : 
VirtualGraphTool12(str,xlabel,ylabel),win(0) {
  WINDOW_TYPE::COLOR_MAP_TYPE *gtkcm=
    dynamic_cast<WINDOW_TYPE::COLOR_MAP_TYPE*>(cm);
  win=OPERATOR_NEW WINDOW_TYPE(str,gtkcm,maxwinsize,maxwinsize,mono);
  if (xxhi<=xxlo) {
    double xlo=min(xxlo,xxhi);
    double xhi=max(xxlo,xxhi);
    double size=exp(floor(log(max(abs(xxlo),abs(xxhi)))));
    xxlo=xlo-size;
    xxhi=xhi+size;
  }
  delayedConstructor(xxlo,xxhi,yylo,yyhi,FALSE,mono);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool::~XYGraphTool() {
  if (win) delete win; win=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XYGraphTool::copyFrom(const XYGraphTool::Buffer &buf) {
  win->refreshFrom(*buf.getWindow());
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XYGraphTool::printOn(ostream &os) const {
  os << "XYGraphTool:" << endl;
  VirtualGraphTool12::printOn(os);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool::Buffer::Buffer(XYGraphTool *gt) :
VirtualGraphTool12(gt->getName(),gt->getXLabel(),gt->getYLabel()),win(0)
{
  WINDOW_TYPE* gtkwin=dynamic_cast<WINDOW_TYPE*>(gt->getWindow());
  win=OPERATOR_NEW WINDOW_TYPE::Buffer(gtkwin);
  delayedConstructor(gt->getLowX(),gt->getHighX(),gt->getLowY(),
    gt->getHighY(),FALSE,gt->monochrome());
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool::Buffer::~Buffer() {
  if (win) delete win; win=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XYGraphTool::Buffer::printOn(ostream &os) const {
  os << "XYGraphTool::Buffer:" << endl;
  VirtualGraphTool12::printOn(os);
}
