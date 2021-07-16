#ifndef _BINOMIAL_COEFFICIENT_H_
#define _BINOMIAL_COEFFICIENT_H_
#include <string.h>
#include "Errors.H"
#include "MemoryDebugger.H"
#include "Tracer.H"

class BinomialCoefficient {
  private:
    int *data;
    int **coefs;
    int order;

    static int numberStoredForOne(int n) { return n/2; }
    static int numberStoredForAll(int n) {
      return (n%2==0 ? n*n : (n-1)*(n+1) )/4;
    }
    explicit BinomialCoefficient(const BinomialCoefficient&);
    BinomialCoefficient& operator=(const BinomialCoefficient&);
  public:
    BinomialCoefficient() : coefs(0),order(-1) {}
    explicit BinomialCoefficient(int);
    ~BinomialCoefficient();
    int value(int n,int k) const {
//    TRACER_CALL(t,"BinomialCoefficient::value");
      CHECK_BOUNDS(n,0,order+1);
      CHECK_BOUNDS(k,0,n+1);
      if (k==0 || k==n) return 1;
      return (k<=n-k ? coefs[n-2][k-1] : coefs[n-2][n-k-1]);
    }
};
#endif
