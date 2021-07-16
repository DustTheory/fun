#include <limits.h>
#include "arch.H"
//#include "LaUtil.H"
int int_zero_=0;
int int_one_=1;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Vector.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#ifdef __GNUC__
  template class Vector<int>;
  template ostream& operator<<(ostream&,const Vector<int> &);
  template void testVector(int);
#endif
template<> int Vector<int>::undefined_ = INT_MAX;
template<> int Vector<int>::zero_ = int_zero_;
template<> int Vector<int>::one_ = int_one_;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Matrix.C"

#ifdef __GNUC__
  template class Matrix<int>;
  template ostream& operator<<(ostream&,const Matrix<int>&);
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Transportation.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#ifdef __GNUC__
  template class Transportation<int>;
  template void testTransportation(int);
#endif
template<> int Transportation<int>::zero_ = 0;
template<> int Transportation<int>::huge_ = INT_MAX;
