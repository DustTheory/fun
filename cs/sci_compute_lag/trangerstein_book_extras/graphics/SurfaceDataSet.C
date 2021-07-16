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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/SurfaceDataSet.C,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#if (SPACEDIM==3)
#include "SurfaceDataSet.H"
//#include "Tracer.H"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//corner numbers:
//    z       3..........7
//    |     . .        . .
//    |   .   .      .   .
//    | .     .    .     .
//    1..........5       .
//    .       .  .       .
//    .       2..........6
//    .     .    .     .
//    .   .      .   .
//    . .        . .
//    0..........4---->x
//edge numbers:
//    z       .....3......
//    |     . .        . .
//    |   6   .      7   .
//    | .     9    .    11
//    ......1.....       .
//    .       .  .       .
//    .       .....2......
//    8     .   10     .
//    .   4      .   5
//    . .        . .
//    ......0.....---->x

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SurfaceDataSet::SurfaceDataSet(VolGraphTool *gt,const Vector3 &l,
const Vector3 &h,int nf,int nc,double (*func)(const Vector3&),double vmin,
double vmax) : graph_tool(gt),low(l),high(h),center((l+h)*HALF),
radius(((h-l)*HALF).max()),var_min(vmin),var_max(vmax),nptsf(nf),
nptsc(nc),surface_func(func) { }
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SurfaceDataSet::plot(bool do_best,const AxisClipPlane *cp) {  
  VolGraphTool::WINDOW_TYPE *glw=
    dynamic_cast<VolGraphTool::WINDOW_TYPE*>(graph_tool->getWindow());
  CHECK_POINTER(glw)

  CHECK_TEST(graph_tool->isDrawing());
  double location=graph_tool->getLocation(cp);
  unsigned int dir=cp->getDirection();
  double fudge=graph_tool->getLength()[dir]*0.0009765625;
  if (cp->getHand()==LEFT_HAND) location += fudge;
  else location -= fudge;
  int d1=(dir+1)%3,d2=(dir+2)%3;

  double dx_val=TWO*radius/double(nptsf);
  Vector3 dx(dx_val,dx_val,dx_val);
  Vector3 centre;
  Vector3 low,high;
  centre[dir]=location; 
  low[dir]=location;
  high[dir]=location;
  low[d1]=center[d1]-radius;
  for (int j=0;j<nptsf;j++) {
    high[d1]=low[d1]+dx[d1];
    centre[d1]=HALF*(low[d1]+high[d1]);
    low[d2]=center[d2]-radius;
    for (int k=0;k<nptsf;k++) {
      high[d2]=low[d2]+dx[d2];
      centre[d2]=HALF*(low[d2]+high[d2]);
      double frac=(surface_func(centre)-var_min)/(var_max-var_min);
      graph_tool->colorRect(frac,low,high,dir);
      low[d2]=high[d2];
    }
    low[d1]=high[d1];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SurfaceDataSet::plot(bool do_best) {  
  VolGraphTool::WINDOW_TYPE *glw=
    dynamic_cast<VolGraphTool::WINDOW_TYPE*>(graph_tool->getWindow());
  CHECK_POINTER(glw)
  VolGraphTool::WINDOW_TYPE::SlideBar *slide_bar=glw->slideBar();
  CHECK_POINTER(slide_bar)
  double frac=graph_tool->isoSurfaceFrac();
  double val=var_min+frac*(var_max-var_min);
  plot(do_best,val);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SurfaceDataSet::plot(bool do_best,double val) {  
  if (var_max<=var_min) return;
  val=max(var_min,min(var_max,val));
  double frac=(val-var_min)/(var_max-var_min);

  int point_number=(do_best ? nptsf : nptsc);
  double dx_val=TWO*radius/double(point_number);
  Vector3 dx(dx_val,dx_val,dx_val);
  Vector3 centre;

  for (int i=0;i<point_number;i++) {
    centre[0]=center[0]+(i+HALF)*dx[0]-radius;
    for (int j=0;j<point_number;j++) {
      centre[1]=center[1]+(j+HALF)*dx[1]-radius;
      for (int k=0;k<point_number;k++) {
        centre[2]=center[2]+(k+HALF)*dx[2]-radius;
        {
          double value[2][2][2];
          Vector3 vertex[2][2][2];
          for (HAND h0=LEFT_HAND;;h0=RIGHT_HAND) {
            for (HAND h1=LEFT_HAND;;h1=RIGHT_HAND) {
              for (HAND h2=LEFT_HAND;;h2=RIGHT_HAND) {
                vertex[h0][h1][h2]=centre+dx*(Vector3(h0,h1,h2)-HALF);
                value[h0][h1][h2]=
                  surface_func(vertex[h0][h1][h2])-val;
                if (h2==RIGHT_HAND) break;
              }
              if (h1==RIGHT_HAND) break;
            }
            if (h0==RIGHT_HAND) break;
          }

          Vector3 intersected[3][2][2];
          for (int dir=0;dir<3;dir++) {
            int id1=(dir+1)%3;
            int id2=(dir+2)%3;
            for (HAND h1=LEFT_HAND;;h1=RIGHT_HAND) {
              for (HAND h2=LEFT_HAND;;h2=RIGHT_HAND) {
                HAND hl[3],hr[3];
                hl[dir]=LEFT_HAND;
                hr[dir]=RIGHT_HAND;
                hl[id1]=hr[id1]=h1;
                hl[id2]=hr[id2]=h2;
                double vl=value[hl[0]][hl[1]][hl[2]];
                double vr=value[hr[0]][hr[1]][hr[2]];
                if (vl*vr<=ZERO) {
                  double alpha=(abs(vl-vr)>ZERO ? vl/(vl-vr) : HALF);
                  intersected[dir][h1][h2]=
                    vertex[hl[0]][hl[1]][hl[2]]*(ONE-alpha)
                   +vertex[hr[0]][hr[1]][hr[2]]*alpha;
                }
                if (h2==RIGHT_HAND) break;
              }
              if (h1==RIGHT_HAND) break;
            }
          }
          graph_tool->drawSurfaceWithinCube(frac,intersected,value);
        }
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SurfaceDataSet::printOn(ostream &os) const {
  os << "SurfaceDataSet: nptsf = " << nptsf << endl;
  os << "\tnptsc = " << nptsc << endl;
  os << "\tlow = " << low << endl;
  os << "\thigh = " << high << endl;
  os << "\tcenter = " << center << endl;
  os << "\tradius = " << radius << endl;
  os << "\tvar_min = " << var_min << endl;
  os << "\tvar_max = " << var_max << endl;
  os << "\tgraph_tool = " << graph_tool << endl;
}
#endif
