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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/VirtualGLWindow.C,v 1.1 2009/08/20 17:31:47 johnt Exp $"
#if (SPACEDIM>1)
#include<string.h>
#include "VirtualGLWindow.H"
//#include "Tracer.H"
# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
MyDisplayList::MyDisplayList() : number(1) {
  base=glGenLists(number);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
MyDisplayList::~MyDisplayList() {
  if (number>0) glDeleteLists(base,number);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void MyDisplayList::call() const {
//see OpenGL Ref Man p 78
  glPushAttrib(GL_ENABLE_BIT);
  if (number==1) glCallList(base); 
  glPopAttrib();
  glFlush();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualGLWindow::VirtualGLWindow(const char *na,const Vector3 &b,
bool rb) : name(0),bound(b),plot_obj(0),rotate_best(rb),
inside_draw_delimiter(false)
#if (SPACEDIM==3)
,number_surfaces(1)
#endif
{
  if (na) {
    name=OPERATOR_NEW_BRACKET(char,strlen(na)+1);
    strcpy(name,na); 
  } else {
    name=OPERATOR_NEW_BRACKET(char,1);
    name[0]='\0';
  }
//clip_planes have to be constructed inside ginit to appear
  for (int i=0;i<3;i++) {
    clip_plane[i][0]=clip_plane[i][1]=0;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualGLWindow::~VirtualGLWindow() { 
  CHECK_TEST(plot_obj==0)
  if (name!=0) delete [] name; name=0;
#ifdef DEBUG
  for (unsigned int i=0;i<3;i++) {
    for (int j=0;j<1;j++) {
      ASSERT(clip_plane[i][j]==0)
    }
  }
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGLWindow::ginit() {
   double fudge=1.-0.03125;
   for (unsigned int i=0;i<3;i++) {
     CHECK_TEST(clip_plane[i][LEFT_HAND]==0);
     CHECK_TEST(clip_plane[i][RIGHT_HAND]==0);
     clip_plane[i][LEFT_HAND]=
       OPERATOR_NEW AxisClipPlane(i,LEFT_HAND,-bound[i]*fudge);
     clip_plane[i][RIGHT_HAND]=
       OPERATOR_NEW AxisClipPlane(i,RIGHT_HAND,bound[i]*fudge);
   }
   enableClipPlanes();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GLdouble VirtualGLWindow::detMatrix(GLdouble *m) const {
  GLdouble det=m[0]*
    (m[5]*m[10]*m[15]
    +m[6]*m[11]*m[13]
    +m[7]*m[ 9]*m[14]
    -m[7]*m[10]*m[13]
    -m[6]*m[ 9]*m[15]
    -m[5]*m[11]*m[14])
              -m[1]*
    (m[4]*m[10]*m[15]
    +m[6]*m[11]*m[12]
    +m[7]*m[ 8]*m[14]
    -m[7]*m[10]*m[12]
    -m[6]*m[ 8]*m[15]
    -m[4]*m[11]*m[14])
              +m[2]*
    (m[4]*m[ 9]*m[15]
    +m[5]*m[11]*m[12]
    +m[7]*m[ 8]*m[13]
    -m[7]*m[ 9]*m[12]
    -m[5]*m[ 8]*m[15]
    -m[4]*m[11]*m[13])
              -m[3]*
    (m[4]*m[ 9]*m[14]
    +m[5]*m[10]*m[12]
    +m[6]*m[ 8]*m[13]
    -m[6]*m[ 9]*m[12]
    -m[5]*m[ 8]*m[14]
    -m[4]*m[10]*m[13]);
  return det;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGLWindow::printMatrix(GLdouble *m) const {
  for (int i=0;i<16;i+=4) {
    cout << "\t" << m[i] << " " << m[i+1] << " " 
         << m[i+2] << " " << m[i+3] << endl;
  }
  cout << "\tdet = " << detMatrix(m) << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGLWindow::clearClipPlanes() { 
  for (unsigned int i=0;i<3;i++) {
    for (int j=0;j<2;j++) {
      if (clip_plane[i][j]!=0) {
        delete clip_plane[i][j]; 
      }
      clip_plane[i][j]=0;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//have to feed user coords to window for cross hairs
void VirtualGLWindow::rescale(const Vector3 &l,const Vector3 &h) {
  CHECK_TEST(l<h);
  user_low=l;
  user_high=h;
  user_center=(user_low+user_high)*HALF;
#if (SPACEDIM==2)
  user_maxlen=max(user_high[0]-user_low[0],user_high[1]-user_low[1]);
#endif
#if (SPACEDIM==3)
  user_maxlen=(user_high-user_low).max();
#endif
  bound=(user_high-user_center)/user_maxlen;
#if (SPACEDIM==2)
  bound[2]=(user_high[2]-user_center[2])/(user_high[2]-user_low[2]);
#endif
  double fudge=1.-0.03125;
  for (int i=0;i<3;i++) {
    double loc=-bound[i]*fudge;
    clip_plane[i][LEFT_HAND]->setLocation(loc);
    loc=bound[i]*fudge;
    clip_plane[i][RIGHT_HAND]->setLocation(loc);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double VirtualGLWindow::clipPlaneLocation(unsigned int d,HAND h) const {
  return clip_plane[d][h]->getLocation();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGLWindow::disableClipPlanes() const {
  for (int i=0;i<3;i++) {
    clip_plane[i][LEFT_HAND]->disable();
    clip_plane[i][RIGHT_HAND]->disable();
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGLWindow::enableClipPlanes() const {
  for (int i=0;i<3;i++) {
    clip_plane[i][LEFT_HAND]->enable();
    clip_plane[i][RIGHT_HAND]->enable();
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGLWindow::debugClipPlanes() const {
  for (int j=0;j<3;j++) {
    cout << "\tclip_plane[" << j <<  "] = "
         << clip_plane[j][LEFT_HAND] << " "
         << clip_plane[j][RIGHT_HAND] << endl;
//  cout << "\tj = " << j << " clip_plane enabled = "
//       << clip_plane[j][LEFT_HAND]->isEnabled() << " "
//       << clip_plane[j][RIGHT_HAND]->isEnabled() << endl;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGLWindow::printOn(ostream &os) const {
//os << "VirtualGLWindow:" 
//   << "\n\tinside_draw_delimiter = " 
//   << inside_draw_delimiter
#if (SPACEDIM==3)
//   << "\n\tnumber_surfaces = " << number_surfaces 
#endif
//   << endl; 
  VirtualWindow::printOn(os);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualGLWindow::ClipPlane::ClipPlane(GLenum pn,GLdouble x,GLdouble y,
GLdouble z,GLdouble w) : plane_number(pn)
#if (SPACEDIM==3)
,display_list(0)
#endif
{
#if (SPACEDIM==3)
  display_list=OPERATOR_NEW MyDisplayList();
#endif
  equation_coefs[0]=x; equation_coefs[1]=y;
  equation_coefs[2]=z; equation_coefs[3]=w;
  fixCoefs();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualGLWindow::DrawDelimiter::DrawDelimiter(bool do_best,
VirtualGLWindow *w) : win(w),is_enabled(true),clipping(false) {
  win->startDrawing(do_best);
}
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualGLWindow::DrawDelimiter::DrawDelimiter(VirtualGLWindow *w,
const ClipPlane *cp) : win(w),clipping(true) {
  if (win!=0) {
    CHECK_TEST(!win->isDrawing())
    win->inside_draw_delimiter=true;
    is_enabled=cp->isEnabled();
    win->disableClipPlanes(); //must be outside glNewList--glEndList
    glNewList(cp->getDisplayList(),GL_COMPILE);
  }
} 
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualGLWindow::DrawDelimiter::~DrawDelimiter() {
  if (win!=0) {
    CHECK_TEST(win->isDrawing())
    if (clipping) {
      win->inside_draw_delimiter=false;
      glEndList();
    } else win->finishDrawing();
    if (is_enabled) win->enableClipPlanes();
  }
}

#include "Pair.H"
#include "Quaternion.H"
template class Pair<double>;
template Pair<double> operator- <> (const Pair<double>&,const Pair<double>&);
template double innerProduct <> (const Pair<double>&,const Pair<double>&);
template Pair<double> closestInteriorPosition(const Pair<double>,
  const Pair<double>,const Pair<double>,const Pair<double>&);
template Quaternion<double> trackBall <> (Pair<double>,Pair<double>);
template ostream& operator<<(ostream&,const Quaternion<double>&);
#endif
