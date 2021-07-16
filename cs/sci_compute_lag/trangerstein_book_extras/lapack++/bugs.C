//#include <complex.h>
#include "f2c.h"
//template <class FLOAT> class complex { };

extern "C" {
  ilaenv_(int *i, const char *n, const char *opts,
    int *n1, int *n2, int *n3, int *n4,
    ftnlen n_len, ftnlen opts_len);
}
