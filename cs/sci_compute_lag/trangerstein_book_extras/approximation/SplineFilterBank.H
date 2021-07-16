#ifndef _SPLINE_FILTER_BANK_H_
#define _SPLINE_FILTER_BANK_H_
#include <string.h>
#include "BinomialCoefficient.H"
#include "Errors.H"
#include "Filter.H"
#include "MemoryDebugger.H"
#include "Tracer.H"

class SplineFilterBank {
  private:
    Filter *dual_highpass;
    Filter *dual_lowpass;
    Filter *highpass;
    Filter *lowpass;
//  prevent public use of copy constructor
    explicit SplineFilterBank(const SplineFilterBank &sfb) { }
    SplineFilterBank& operator=(const SplineFilterBank &sfb) {
      return *this;
    }
  public:
    SplineFilterBank() : dual_highpass(0),dual_lowpass(0),highpass(0),
      lowpass(0) {}
    explicit SplineFilterBank(int,int);
    ~SplineFilterBank() {
      delete lowpass;
      delete highpass;
      delete dual_lowpass;
      delete dual_highpass;
    }
    inline const Filter& dualHighpassFilter() const {
      return *dual_highpass;
    }
    inline const Filter& dualLowpassFilter() const { return *dual_lowpass; }
    inline const Filter& highpassFilter() const { return *highpass; }
    inline const Filter& lowpassFilter() const { return *lowpass; }
    int analysisFilterStart() const {
      return min(dual_lowpass->firstIndex(),dual_highpass->firstIndex());
    }
    int analysisFilterFinish() const {
      return max(dual_lowpass->lastIndex(),dual_highpass->lastIndex());
    }
    int synthesisFilterStart() const {
      return min(lowpass->firstIndex(),highpass->firstIndex());
    }
    int synthesisFilterFinish() const {
      return max(lowpass->lastIndex(),highpass->lastIndex());
    }
    void analyze(const Signal&,Signal&,Signal&) const;
    void synthesize(const Signal&,const Signal&,Signal&) const;
};
#endif
