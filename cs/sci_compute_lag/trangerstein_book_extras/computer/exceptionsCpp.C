#include <cmath> // for isinf, isfinite
#include <cstdlib> // for EXIT_SUCCESS
#include <iostream> // for cout
#include <iomanip> // for setw, setprecision
#include <limits> // for signaling_NaN, quiet_NaN, infinity
#define __USE_GNU 1
#include <fenv.h>

#if (_TRAP_FPE_==1)
static void __attribute__((constructor)) trapfpe () {
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}
#endif

using namespace std;
//from bits/fenv.h:
//  FE_INVALID = 0x1
//  __FE_DENORM = 0x2
//  FE_DIVBYZERO = 0x4
//  FE_OVERFLOW = 0x8
//  FE_UNDERFLOW = 0x10
//  FE_INEXACT = 0x20

int main(int /* argc */,char** /* argv */) {
  cout << boolalpha;
  cout << setprecision(17) << hex;
  cout << "exception codes:\n\t invalid = " << FE_INVALID 
    << "\n\t denorm = = " << __FE_DENORM << "\n\t divide by zero = "
    << FE_DIVBYZERO << endl;
  cout << "\t overflow = " << FE_OVERFLOW << "\n\t underflow = "
    << FE_UNDERFLOW << "\n\t inexact = " << FE_INEXACT << endl;

  try { // no exception
    double d=numeric_limits<double>::infinity();
    cout << "numeric_limits<double>::infinity() = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "d=numeric_limits<double>::infinity() threw exception "
         << e.what() << endl;
  } catch (...) {
    cout << "d=numeric_limits<double>::infinity() threw unknown exception"
         << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

#if (_TRAP_FPE_!=1)
  try { // FE_INVALID generated but no exception thrown
        // floating point exception will be trapped
    double d=std::numeric_limits<double>::signaling_NaN();
    cout << "signaling_NaN = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "d=signaling_NaN() threw exception " << e.what() << endl;
  } catch (...) {
    cout << "d=signaling_NaN() threw unknown exception" << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  try { // no exception
    double d=std::numeric_limits<double>::quiet_NaN();
    cout << "quiet_NaN = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "d=quiet_NaN() threw exception " << e.what() << endl;
  } catch (...) {
    cout << "d=quiet_NaN() threw unknown exception" << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

#if (_TRAP_FPE_!=1)
  try { // no exception found in test
        // floating point exception will be trapped
    double d=1.;
    d /= numeric_limits<double>::infinity();
    cout << "1./numeric_limits<double>::infinity() = " << setw(24) << d
         << " " << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | __FE_DENORM | FE_DIVBYZERO
           | FE_OVERFLOW | FE_UNDERFLOW | FE_INEXACT) << endl;
  } catch (exception &e) {
    cout << "1./numeric_limits<double>::infinity() threw exception "
         << e.what() << endl;
  } catch (...) {
    cout
      << "1./numeric_limits<double>::infinity()() threw unknown exception"
      << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  try { // FE_OVERFLOW generated but no exception thrown
        // floating point exception will be trapped
    double d=numeric_limits<double>::max();
    d *= 2.;
    cout << "numeric_limits<double>::max() * 2. = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "numeric_limits<double>::max() * 2. threw exception "
         << e.what() << endl;
  } catch (...) {
    cout << "numeric_limits<double>::max() * 2. threw unknown exception "
         << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  try { // no exception
    double d=numeric_limits<double>::infinity();
    d *= numeric_limits<double>::infinity();
    cout << "infinity * infinity = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "infinity * infinity threw exception " << e.what()
         << endl;
  } catch (...) {
    cout << "infinity * infinity threw unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

#if (_TRAP_FPE_!=1)
  try { // FE_DIVBYZERO generated but no exception thrown
        // floating point exception will be trapped
    double d=1.;
    double z=0.;
    d /= z;
    cout << "1. / 0. = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "1. / 0. threw exception " << e.what() << endl;
  } catch (...) {
    cout << "1. / 0. threw non-standard exception" << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  try { // FE_INVALID generated but no exception thrown
        // floating point exception will be trapped
    double d=0.;
    double z=0.;
    d /= z;
    cout << "0. / 0. = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "0. / 0. threw exception " << e.what() << endl;
  } catch (...) {
    cout << "0. / 0. threw unknown unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  try { // no exception
    double d=numeric_limits<double>::infinity();
    d += numeric_limits<double>::infinity();
    cout << "infinity + infinity = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "infinity + infinity threw exception " << e.what() << endl;
  } catch (...) {
    cout << "infinity + infinity threw unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

#if (_TRAP_FPE_!=1)
  try { // FE_INVALID generated but no exception thrown
        // floating point exception will be trapped
    double d=numeric_limits<double>::infinity();
    d -= numeric_limits<double>::infinity();
    cout << "infinity - infinity = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "infinity - infinity threw exception " << e.what() << endl;
  } catch (...) {
    cout << "infinity - infinity threw unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  try { // FE_INVALID generated but no exception thrown
        // floating point exception will be trapped
    double d=numeric_limits<double>::infinity();
    d /= numeric_limits<double>::infinity();
    cout << "infinity / infinity = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "infinity / infinity threw exception " << e.what() << endl;
  } catch (...) {
    cout << "infinity / infinity threw unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  try { // FE_INVALID generated but no exception thrown
        // floating point exception will be trapped
    double d=numeric_limits<double>::infinity();
    d *= 0.;
    cout << "infinity * 0. = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "infinity * 0. threw exception " << e.what() << endl;
  } catch (...) {
    cout << "infinity * 0. threw unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  try { // FE_INVALID generated but no exception thrown
        // floating point exception will be trapped
    double d=log(-1.);
    cout << "log(-1.) = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "log(-1.) threw exception " << e.what() << endl;
  } catch (...) {
    cout << "log(-1.) threw unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  try { // FE_DIVBYZERO generated but no exception thrown
        // floating point exception will be trapped
    double d=log(0.);
    cout << "log(0.) = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "log(0.) threw exception " << e.what() << endl;
  } catch (...) {
    cout << "log(0.) threw unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  try { // FE_INVALID generated but no exception thrown
        // floating point exception will be trapped
    double d=sqrt(-1.);
    cout << "sqrt(-1.) = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "sqrt(-1.) threw exception " << e.what() << endl;
  } catch (...) {
    cout << "sqrt(-1.) threw unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  try { // FE_INVALID generated but no exception thrown
        // floating point exception will be trapped
    double d=acos(2.);
    cout << "acos(2.) = " << setw(24) << d << " "
         << *reinterpret_cast<long long*>(&d) << " "
         << fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) << endl;
  } catch (exception &e) {
    cout << "acos(2.) threw exception " << e.what() << endl;
  } catch (...) {
    cout << "acos(2.) threw unknown exception " << endl;
  }
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  return EXIT_SUCCESS;
}
