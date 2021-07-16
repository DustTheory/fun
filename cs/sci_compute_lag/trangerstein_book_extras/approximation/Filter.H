#ifndef _FILTER_H_
#define _FILTER_H_
#include <string.h>
#include "Errors.H"
#include "MemoryDebugger.H"
#include "Signal.H"
#include "Tracer.H"

class Filter {
  private:
    Signal *impulse_response;

//  prevent public use of copy
    Filter& operator=(const Filter &filter) {
//    impulse_response=filter.impulse_response;
//    return *this;
    }
  public:
    Filter() : impulse_response(0) {}
    explicit Filter(int f,int l) {
      impulse_response=OPERATOR_NEW Signal(f,l);
    }
    explicit Filter(const Filter &f) {
      impulse_response=OPERATOR_NEW Signal(*f.impulse_response);
    }
    ~Filter() {
      delete impulse_response; impulse_response=0;
    }
    inline int firstIndex() const { return impulse_response->firstIndex(); }
    inline int lastIndex() const { return impulse_response->lastIndex(); }
    inline double value(int n) const {
      return impulse_response->value(n);
    }
    inline double& value(int n) {
      return impulse_response->value(n);
    }
    inline const Signal& impulseResponse() const { 
      return *impulse_response;
    }
    inline Signal& impulseResponse() { return *impulse_response; }
    inline void replaceWithAdjoint() { 
      impulse_response->replaceWithAdjoint();
    }
    void operateOn(const Signal&,Signal&) const;
    void operateAdjointOn(const Signal&,Signal&) const;
    void operateAdjointOnAndDownsample(const Signal&,Signal&) const;
    void addUpsampleAndOperateOn(const Signal&,Signal&) const;
};
#endif
