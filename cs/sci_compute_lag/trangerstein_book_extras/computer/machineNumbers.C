//#include <cmath> // for isinf, isfinite
#include <cstdlib> // for EXIT_SUCCESS
#include <iostream> // for cout
#include <iomanip> // for setw, setprecision
#include <limits> // for numeric_limits
using namespace std;

int main(int /* argc */,char** /* argv */) {
  cout << boolalpha;
  cout << setprecision(17) << hex;

  float f=numeric_limits<float>::epsilon();
  cout << "\tnumeric_limits<float>::epsilon() = " << setw(16) << f << " "
       << *reinterpret_cast<int*>(&f) << endl;
  f=numeric_limits<float>::min();
  cout << "\tnumeric_limits<float>::min() = " << setw(16) << f << " "
       << *reinterpret_cast<int*>(&f) << endl;
  f=numeric_limits<float>::max();
  cout << "\tnumeric_limits<float>::max() = " << setw(16) << f << " "
       << *reinterpret_cast<int*>(&f) << endl;
  f=numeric_limits<float>::infinity();
  cout << "\tnumeric_limits<float>::infinity() = " << setw(16) << f << " "
       << *reinterpret_cast<int*>(&f) << endl;
  cout << endl;

  double d=numeric_limits<double>::epsilon();
  cout << "numeric_limits<double>::epsilon() = " << setw(24) << d << " "
       << *reinterpret_cast<long long*>(&d) << endl;
  d=numeric_limits<double>::min();
  cout << "numeric_limits<double>::min()= " << setw(24) << d << " "
       << *reinterpret_cast<long long*>(&d) << endl;
  d=numeric_limits<double>::max();
  cout << "numeric_limits<double>::max() = " << setw(24) << d << " "
       << *reinterpret_cast<long long*>(&d) << endl;
  d=numeric_limits<double>::infinity();
  cout << "numeric_limits<double>::infinity() = " << setw(24) << d << " "
       << *reinterpret_cast<long long*>(&d) << endl;
  cout << endl;

  d=M_E;
  cout << "M_E = " << setw(24) << d << endl;
  d=M_LOG10E;
  cout << "M_LOG10E = " << setw(24) << d << endl;
  d=M_LN2;
  cout << "M_LN2 = " << setw(24) << d << endl;
  d=M_PI;
  cout << "M_PI = " << setw(24) << d << endl;
  d=M_PI_2;
  cout << "M_PI_2 = " << setw(24) << d << endl;
  d=M_PI_4;
  cout << "M_PI_4 = " << setw(24) << d << endl;
  d=M_1_PI;
  cout << "M_1_PI = " << setw(24) << d << endl;
  d=M_2_PI;
  cout << "M_2_PI = " << setw(24) << d << endl;
  d=M_SQRT2;
  cout << "M_SQRT2 = " << setw(24) << d << endl;
  d=M_SQRT1_2;
  cout << "M_SQRT1_2 = " << setw(24) << d << endl;
  cout << endl;

  short s=SHRT_MIN;
  cout << "SHRT_MIN = " << setw(20) << dec << s << " " << hex
       << *reinterpret_cast<short*>(&s) << endl;
  s=SHRT_MAX;
  cout << "SHRT_MAX = " << setw(20) << dec << s << " " << hex
       << *reinterpret_cast<short*>(&s) << endl;

  int i=INT_MIN;
  cout << "INT_MIN = " << setw(20) << dec << i << " " << hex
       << *reinterpret_cast<int*>(&i) << endl;
  i=INT_MAX;
  cout << "INT_MAX = " << setw(20) << dec << i << " " << hex
       << *reinterpret_cast<int*>(&i) << endl;

  long l=LONG_MIN;
  cout << "LONG_MIN = " << setw(20) << dec << l << " " << hex
       << *reinterpret_cast<long*>(&l) << endl;
  l=LONG_MAX;
  cout << "LONG_MAX = " << setw(20) << dec << l << " " << hex
       << *reinterpret_cast<long*>(&l) << endl;

  long long ll=LLONG_MIN;
  cout << "LLONG_MIN = " << setw(20) << dec << ll << " " << hex
       << *reinterpret_cast<long long*>(&ll) << endl;
  ll=LLONG_MAX;
  cout << "LLONG_MAX = " << setw(20) << dec << ll << " " << hex
       << *reinterpret_cast<long long*>(&ll) << endl;
  return EXIT_SUCCESS;
}
