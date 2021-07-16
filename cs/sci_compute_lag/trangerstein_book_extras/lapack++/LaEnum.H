#ifndef LA_ENUM_H
#define LA_ENUM_H "%W%\t%G%"
#include <iostream>
using namespace std;
//until we can figure out how to make these into class members,
//  we are making them global
enum TRANSPOSE_OPTION { NO_MATRIX_TRANSPOSE, MATRIX_TRANSPOSE };
enum PIVOT_OPTION {NO_PIVOTING,PIVOT_ROWS,PIVOT_ROWS_AND_COLUMNS };
enum STATUS_OPTION{ UNKNOWN, INFEASIBLE, PRIMAL_FEASIBLE, DUAL_FEASIBLE,
                    UNBOUNDED, OPTIMAL};

ostream& operator<<(ostream&,PIVOT_OPTION);
ostream& operator<<(ostream&,TRANSPOSE_OPTION);
ostream& operator<<(ostream&,STATUS_OPTION);
#endif
