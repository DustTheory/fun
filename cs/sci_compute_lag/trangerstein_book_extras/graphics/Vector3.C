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
// "$Header: /home/faculty/johnt/cvs/graphics/Vector3.C,v 1.2 2007/01/10 13:18:54 johnt Exp $"
#include "Vector3.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ostream& operator<<(ostream &os, const Vector3 &p) {
  os << '<' << p[0] << ',' << p[1] << ',' << p[2] << '>';
  return os;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ofstream& operator<<(ofstream &os, const Vector3 &p) {
  REAL v[3]; p.getVect(v);
  os.write(reinterpret_cast<char*>(v), SPACEDIM*sizeof(REAL));
  return os;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ifstream& operator>>(ifstream &is, Vector3 &p) {
  is.read(reinterpret_cast<char*>(p.vect), 3*sizeof(REAL));
  return is;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#define STRIP(c) while(is.get() != c)
istream& operator>>(istream &is, Vector3 &p) {
  STRIP('<'); 
  is >> p[0]; STRIP(','); is >> p[1]; STRIP(','); is >> p[2]; 
  STRIP('>');
  return is;
}
#undef STRIP
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Vector3::printOn( ostream& os ) const {
  os << "Vecto3r: " << *this << '\n';
}
