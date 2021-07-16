#ifndef _DAUBECHIES_FILTER_BANK_H_
#define _DAUBECHIES_FILTER_BANK_H_
#include <string.h>
#include "Filter.H"

class DaubechiesFilterBank {
  private:
    Filter *highpass;
    Filter *lowpass;
//  prevent public use of copy constructor
    explicit DaubechiesFilterBank(const DaubechiesFilterBank &sfb) { }
    DaubechiesFilterBank& operator=(const DaubechiesFilterBank &sfb) {
      return *this;
    }
  public:
    DaubechiesFilterBank() : highpass(0),lowpass(0) {}
    explicit DaubechiesFilterBank(int);
    ~DaubechiesFilterBank() {
      delete lowpass;
      delete highpass;
    }
    inline const Filter& dualHighpassFilter() const {
      return *highpass;
    }
    inline const Filter& dualLowpassFilter() const { return *lowpass; }
    inline const Filter& highpassFilter() const { return *highpass; }
    inline const Filter& lowpassFilter() const { return *lowpass; }
    int analysisFilterStart() const {
      return min(lowpass->firstIndex(),highpass->firstIndex());
    }
    int analysisFilterFinish() const {
      return max(lowpass->lastIndex(),highpass->lastIndex());
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
