#include <cmath> // for isinf, isfinite
#include <cstdlib> // for EXIT_SUCCESS
#include <iostream> // for cout
#include <iomanip> // for setw, setprecision
#include <limits> // for numeric_limits
#include <stdio.h> // for printf
using namespace std;

//since no integer type has as many bytes as long double, have to cheat:
typedef union {
  long double ld;
  char bytes[sizeof(long double)];
} converter; 
void print_converter(converter c) {
  printf("0x");
  for (int i=0;i<sizeof(long double);i++) printf("%x",c.bytes[i]);
}

int main(int /* argc */,char** /* argv */) {
  cout << boolalpha;
  cout << "\tsizeof( float ) = " << sizeof(float) << endl;
  cout << "\tsizeof( double ) = " << sizeof(double) << endl;
  cout << "\tsizeof( long double ) = " << sizeof(long double) << endl;
#if (__BYTE_ORDER == __LITTLE_ENDIAN)
  cout << "__BYTE_ORDER == __LITTLE_ENDIAN " << endl;
#else
  cout << "__BYTE_ORDER == __BIG_ENDIAN = " << endl;
#endif
  cout << setprecision(24) << hex;

  float f=numeric_limits<float>::infinity();
  cout << "numeric_limits<float>::infinity() = " << setw(32) << f << " "
       << *reinterpret_cast<int*>(&f) << endl;
  f=numeric_limits<float>::max();
  cout << "numeric_limits<float>::max() = " << setw(32) << f << " "
       << *reinterpret_cast<int*>(&f) << endl;
  f=numeric_limits<float>::min();
  cout << "numeric_limits<float>::min() = " << setw(32) << f << " "
       << *reinterpret_cast<int*>(&f) << endl;
  f=numeric_limits<float>::epsilon();
  cout << "numeric_limits<float>::epsilon() = " << setw(32) << f << " "
       << *reinterpret_cast<int*>(&f) << endl;
  cout << endl;

  cout << setprecision(48) << hex;
  double d=numeric_limits<double>::infinity();
  cout << "numeric_limits<double>::infinity() = " << setw(56) << d << " "
       << *reinterpret_cast<long*>(&d) << endl;
  d=numeric_limits<double>::max();
  cout << "numeric_limits<double>::max() = " << setw(56) << d << " "
       << *reinterpret_cast<long*>(&d) << endl;
  d=numeric_limits<double>::min();
  cout << "numeric_limits<double>::min() = " << setw(56) << d << " "
       << *reinterpret_cast<long*>(&d) << endl;
  d=numeric_limits<double>::epsilon();
  cout << "numeric_limits<double>::epsilon() = " << setw(56) << d << " "
       << *reinterpret_cast<long*>(&d) << endl;
  cout << endl;

  converter c;
  cout << setprecision(112) << hex;
  c.ld=numeric_limits<long double>::infinity();
  cout << "numeric_limits<double>::infinity() = " << setw(120) << c.ld
       << "  ";
  print_converter(c);
  cout << endl;
  c.ld=numeric_limits<long double>::max();
  cout << "numeric_limits<double>::max() = " << setw(120) << c.ld << "  ";
  print_converter(c);
  cout << endl;
  c.ld=numeric_limits<long double>::min();
  cout << "numeric_limits<double>::min() = " << setw(120) << c.ld << "  ";
  print_converter(c);
  cout << endl;
  c.ld=numeric_limits<long double>::epsilon();
  cout << "numeric_limits<double>::epsilon() = " << setw(120) << c.ld
       << "  ";
  print_converter(c);
  cout << endl;

//rounding error occurs after mantissa bits become all 1's
  cout << setprecision(24) << hex;
  cout << "compute a_n = 0.5*a_{n-1} + 1., a_1 = 1 until rounding error"
       << endl;
  float aold=1.;
  float anew=1.+0.5*aold;
  while (anew < 2.) {
    aold=anew;
    cout << setw(24) << anew << " " << *reinterpret_cast<int*>(&anew)
         << endl;
    anew=1.+0.5*aold;
  }
  cout << setw(24) << anew << " " << *reinterpret_cast<int*>(&anew)
       << endl;

//overflow
//bits in mantissa ordered most significant first
//rounding error occurs when number becomes even
//  the rounding error turns number into a power of 2
  cout << endl;
  cout << "compute a_n = 2.*a_{n-1} + 1., a_1 = 1 until overflow" << endl;
  float a=1.;
  while (isfinite(a)) {
    cout << setw(12) << a << " " << *reinterpret_cast<int*>(&a) << endl;
    a=2.*a+1.;
  }
  cout << setw(12) << a << " " << *reinterpret_cast<int*>(&a) << endl;

//underflow
//de-normal numbers occur when biased exponent becomes zero, and
//  then the mantissa bits are no longer all ones
  cout << endl;
  cout << "compute a_n = 0.5*a_{n-1}, a_1 = 2-epsilon while positive"
       << endl;
  while (aold > 0.) {
    cout << setw(12) << aold << " " << *reinterpret_cast<int*>(&aold)
         << endl;
    aold=0.5*aold;
  }
  cout << setw(12) << aold << " " << *reinterpret_cast<int*>(&aold) << endl;
  return EXIT_SUCCESS;
}
