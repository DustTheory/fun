#ifndef _SIGNAL_H_
#define _SIGNAL_H_
#include "Errors.H"

class Signal {
  private:
    double *data;
    int first,last;
  public:
    Signal() : data(0),first(0),last(-1) {}
    explicit Signal(int,int);
    explicit Signal(const Signal&);
    ~Signal();
    Signal& operator=(double);
    Signal& operator*=(double);
    inline int firstIndex() const { return first; }
    inline int lastIndex() const { return last; }
    inline int length() const { return last-first+1; }
    inline double value(int n) const { 
      CHECK_BOUNDS(n,first,last+1)
      return data[n-first];
    }
    inline double& value(int n) { 
      CHECK_BOUNDS(n,first,last+1)
      return data[n-first];
    }
    inline void shift(int s) { first-=s; last-=s; }
    void reverse();
    void replaceWithAdjoint();
};
#endif
