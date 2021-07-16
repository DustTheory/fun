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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/AreaGraphTool.C,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#if (SPACEDIM>1)
#include "AreaGraphTool.H"
#include <algorithm>
//#include "Tracer.H"
# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool::AreaGraphTool(const char *str,const char *xl,
const char *yl,const double *rlo,const double *rhi,VirtualColormap *cm,
char *dn,double maxwinsize,bool mono) : VirtualGraphTool12(str,xl,yl),
win(0) {
  double maxlen=max(rhi[0]-rlo[0],rhi[1]-rlo[1]);
  double width=maxwinsize*(rhi[0]-rlo[0])/maxlen;
  double height=maxwinsize*(rhi[1]-rlo[1])/maxlen;
  WINDOW_TYPE::COLOR_MAP_TYPE *gtkcm=
    dynamic_cast<WINDOW_TYPE::COLOR_MAP_TYPE*>(cm);
  win=OPERATOR_NEW WINDOW_TYPE(str,gtkcm,width,height);
  delayedConstructor(rlo[0],rhi[0],rlo[1],rhi[1],true,mono);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool::AreaGraphTool(const char *str,const char *xl,
const char *yl,VirtualColormap *cm,char *dn,double width,double height,
bool mono) : VirtualGraphTool12(str,xl,yl),win(0) {
  WINDOW_TYPE::COLOR_MAP_TYPE *gtkcm=
    dynamic_cast<WINDOW_TYPE::COLOR_MAP_TYPE*>(cm);
  win=OPERATOR_NEW WINDOW_TYPE(str,gtkcm,width,height);
  delayedConstructor(width,height,mono);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void AreaGraphTool::printOn(ostream& os) const {
  os << "AreaGraphTool: win = " << win << endl;
  win->printOn(os);
  VirtualGraphTool12::printOn(os);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool::Buffer::Buffer(AreaGraphTool *gt) :
VirtualGraphTool12(gt->getName(),gt->getXLabel(),gt->getYLabel()),win(0)
{
  WINDOW_TYPE* gtkwin=dynamic_cast<WINDOW_TYPE*>(gt->getWindow());
  win=OPERATOR_NEW WINDOW_TYPE::Buffer(gtkwin);
  delayedConstructor(gt->getLowX(),gt->getHighX(),gt->getLowY(),
    gt->getHighY(),false,gt->monochrome());
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void AreaGraphTool::Buffer::colorQuad(double ratio,double *xpos,double *ypos){
  double xfrac[4],yfrac[4];
  for (int i=0;i<4;i++) {
    xfrac[i]=(xpos[i]-xlo)/(xhi-xlo);
    yfrac[i]=(ypos[i]-ylo)/(yhi-ylo);
  }
  int four=4;
  getWindow()->colorPolygon(&ratio,&four,xfrac,yfrac);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void AreaGraphTool::Buffer::printOn(ostream &os) const {
  os << "AreaGraphTool::Buffer: parent = " << parent << endl;
  parent->printOn(os);
  os << "\n\twin = " << win << endl;
  win->printOn(os);
  VirtualGraphTool12::printOn(os);
}
#endif
