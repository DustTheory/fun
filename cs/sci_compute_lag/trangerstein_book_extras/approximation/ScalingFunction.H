#ifndef _SCALING_FUNCTION_H_
#define _SCALING_FUNCTION_H_
#include "Signal.H"
//#include "Types.H"

class ScalingFunction {
  private:
    const Signal *lowpass;
    Signal *values_at_integers;

    explicit ScalingFunction(const ScalingFunction&);
  public:
    ScalingFunction() : lowpass(0),values_at_integers(0) {}
    explicit ScalingFunction(const Signal&);
    ~ScalingFunction();
    inline double supportLow() const {
      CHECK_POINTER(lowpass);
      return static_cast<double>(lowpass->firstIndex());
    }
    inline double supportHigh() const {
      CHECK_POINTER(lowpass);
      return static_cast<double>(lowpass->lastIndex());
    }
    inline const Signal& signal() const { 
      CHECK_POINTER(lowpass);
      return *lowpass;
    }
    Signal* values(int) const;
};
#endif
