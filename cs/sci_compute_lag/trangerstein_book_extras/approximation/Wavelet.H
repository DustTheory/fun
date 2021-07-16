#ifndef _WAVELET_H_
#define _WAVELET_H_
//#include "NumPtr.H"
#include "ScalingFunction.H"
//#include "Types.H"
//#include "Vector.H"

class Wavelet {
  private:
    const Signal *highpass;
    const ScalingFunction *sf;

    explicit Wavelet(const Wavelet&) : highpass(0),sf(0) {}
  public:
    Wavelet() : highpass(0),sf(0) {}
    explicit Wavelet(const Signal &sig,const ScalingFunction &scfun) :
      highpass(&sig),sf(&scfun) {}
    ~Wavelet() {}
    inline double supportLow() const {
      CHECK_POINTER(highpass);
      return 0.5*(static_cast<double>(highpass->firstIndex())
                 +sf->supportLow());
    }
    inline double supportHigh() const {
      CHECK_POINTER(highpass);
      return 0.5*(static_cast<double>(highpass->lastIndex())
                 +sf->supportHigh());
    }
    inline const Signal& signal() const {
      CHECK_POINTER(highpass);
      return *highpass;
    }
    Signal* values(int) const;
};
#endif
