#include "LaEnum.H"

ostream& operator<<(ostream &s,TRANSPOSE_OPTION tr) {
  switch (tr) {
    case NO_MATRIX_TRANSPOSE:
      s << "no_matrix_transpose";
      break;
    default:
      s << "matrix_transpose";
  }
  return s;
}

ostream& operator<<(ostream &s,PIVOT_OPTION po) {
  switch (po) {
    case NO_PIVOTING:
      s << "no_pivoting";
      break;
    case PIVOT_ROWS_AND_COLUMNS:
      s << "pivot_rows_and_columns";
      break;
    default:
      s << "pivot_rows";
  }
  return s;
}

ostream& operator<<(ostream &s,STATUS_OPTION st) {
  switch (st) {
    case INFEASIBLE:
      s << "infeasible";
      break;
    case PRIMAL_FEASIBLE:
      s << "primal_feasible";
      break;
    case DUAL_FEASIBLE:
      s << "dual_feasible";
      break;
    case UNBOUNDED:
      s << "unbounded";
      break;
    case OPTIMAL:
      s << "optimal";
      break;
    default:
      s << "unknown";
  }
  return s;
}
