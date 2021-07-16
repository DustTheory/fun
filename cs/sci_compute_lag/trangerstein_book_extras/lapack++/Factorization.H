#ifndef _FACTORIZATION_H
#define _FACTORIZATION_H
#include <iostream>
using namespace std;
//#include "SquareMatrix.H"
//#include "TrapezoidalMatrix.H"

namespace Factorization {
  enum SIDE_OPTION {LEFT_SIDE,RIGHT_SIDE};
  enum TRANSPOSE_OPTION {NO_TRANSPOSE,TRANSPOSE,CONJUGATE_TRANSPOSE};
  enum PIVOT_OPTION {NO_PIVOTING,PIVOT_ROWS,PIVOT_ROWS_AND_COLUMNS,
    PIVOT_COLUMNS};
  enum CONDITION_NUMBER_NORM {ONE_NORM,INFINITY_NORM};
  enum EQUILIBRATE_OPTION{NO_EQUILIBRATION,EQUILIBRATE_ROWS,
    EQUILIBRATE_COLUMNS,EQUILIBRATE_ROWS_AND_COLUMNS}; // dlaqge
    // general: dgebal,dgeequ, dgeequb, dlaqge
    // band: dgbequb, dgbequ, dlaqgb
};
//ostream& operator<<(ostream&,Factorization::SIDE_OPTION);
//ostream& operator<<(ostream&,Factorization::TRANSPOSE_OPTION);
//ostream& operator<<(ostream&,Factorization::PIVOT_OPTION);
//ostream& operator<<(ostream&,Factorization::EQUILIBRATE_OPTION);
//Factorization::EQUILIBRATE_OPTION
//  operator++(Factorization::EQUILIBRATE_OPTION&,int);
//Factorization::PIVOT_OPTION
//  operator++(Factorization::PIVOT_OPTION&,int);

inline ostream& operator<<(ostream &os,Factorization::SIDE_OPTION so) {
  os << (so==Factorization::LEFT_SIDE ? "LEFT_SIDE" : "RIGHT_SIDE");
  return os;
}

inline ostream& operator<<(ostream &os,Factorization::TRANSPOSE_OPTION to) {
  switch (to) {
    case Factorization::TRANSPOSE:
      os << "TRANSPOSE";
      break;
    case Factorization::CONJUGATE_TRANSPOSE:
      os << "CONJUGATE_TRANSPOSE";
      break;
    default:
      os << "NO_TRANSPOSE";
  }
  return os;
}

inline ostream& operator<<(ostream &os,Factorization::PIVOT_OPTION po) {
  switch (po) {
    case Factorization::PIVOT_ROWS:
      os << "PIVOT_ROWS";
      break;
    case Factorization::PIVOT_ROWS_AND_COLUMNS:
      os << "PIVOT_ROWS_AND_COLUMNS";
      break;
    case Factorization::PIVOT_COLUMNS:
      os << "PIVOT_COLUMNS";
      break;
    default:
      os << "NO_PIVOTING";
  }
  return os;
}

inline ostream& operator<<(ostream &os,Factorization::EQUILIBRATE_OPTION eo)
{
  switch (eo) {
    case Factorization::EQUILIBRATE_ROWS:
      os << "EQUILIBRATE_ROWS";
      break;
    case Factorization::EQUILIBRATE_COLUMNS:
      os << "EQUILIBRATE_COLUMNS";
      break;
    case Factorization::EQUILIBRATE_ROWS_AND_COLUMNS:
      os << "EQUILIBRATE_ROWS_AND_COLUMNS";
      break;
    default:
      os << "NO_EQUILIBRATION";
  }
  return os;
}

inline ostream& operator<<(ostream &os,
Factorization::CONDITION_NUMBER_NORM cnn) {
  switch (cnn) {
    case Factorization::ONE_NORM:
      os << "ONE_NORM";
      break;
    case Factorization::INFINITY_NORM:
      os << "INFINITY_NORM";
      break;
  }
  return os;
}

inline Factorization::EQUILIBRATE_OPTION
operator++(Factorization::EQUILIBRATE_OPTION &eo,int) {
  switch (eo) {
    case Factorization::NO_EQUILIBRATION:
      return eo=Factorization::EQUILIBRATE_ROWS;
    case Factorization::EQUILIBRATE_ROWS:
      return eo=Factorization::EQUILIBRATE_COLUMNS;
    default:
      return eo=Factorization::EQUILIBRATE_ROWS_AND_COLUMNS;
  }
}

inline Factorization::PIVOT_OPTION
operator++(Factorization::PIVOT_OPTION &po,int) {
  return (po==Factorization::NO_PIVOTING ? po=Factorization::PIVOT_ROWS
    : po=Factorization::PIVOT_ROWS_AND_COLUMNS);
}
#endif 
