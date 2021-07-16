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
// "$Header: /home/faculty/johnt/cvs/graphics/DVector.C,v 1.2 2007/01/10 13:18:54 johnt Exp $"
#include "DVector.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ostream& operator<<(ostream &os, const DVector &p) {
#if (SPACEDIM==1)
  os << '<' << p[0] << '>';
#endif
#if (SPACEDIM==2)
  os << '<' << p[0] << ',' << p[1] << '>';
#endif
#if (SPACEDIM==3)
  os << '<' << p[0] << ',' << p[1] << ',' << p[2] << '>';
#endif
  return os;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ofstream& operator<<(ofstream &os, const DVector &p) {
  REAL v[SPACEDIM]; p.getVect(v);
  os.write(reinterpret_cast<char*>(v), SPACEDIM*sizeof(REAL));
  return os;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ifstream& operator>>(ifstream &is, DVector &p) {
  is.read(reinterpret_cast<char*>(p.vect), SPACEDIM*sizeof(REAL));
  return is;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#define STRIP(c) while(is.get() != c)
istream& operator>>(istream &is, DVector &p) {
  STRIP('<'); 
#if (SPACEDIM==1)
  is >> p[0]; 
#endif
#if (SPACEDIM==2)
  is >> p[0]; STRIP(','); is >> p[1]; 
#endif
#if (SPACEDIM==3)
  is >> p[0]; STRIP(','); is >> p[1]; STRIP(','); is >> p[2]; 
#endif
  STRIP('>');
  return is;
}
#undef STRIP
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DVector::printOn( ostream& os ) const {
  os << "DVector: " << *this << '\n';
}
