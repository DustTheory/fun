#include <cmath> // cos, sin
#include <cstdlib> // for EXIT_SUCCESS
#include <fstream> // for ofstream
#include <iostream> // for cout
#include <iomanip> // for setw, setprecision
#include <limits> // for numeric_limits

using namespace std;

inline double f(double x) { return sin(x); }
inline double fprime(double x) { return cos(x); }

int main(int argc,char** argv) {
  cout << "enter x" << endl;
  double x;
  cin >> x;
  double h=1.;
  int nlevels=64;

  double exact=fprime(x);
  ofstream out_file;
  out_file.open("cpp_derivative_output",ios::out);
  for (int i=0;i<nlevels;i++,h*=0.5) {
    double diff=( f(x+h) - f(x) )/h;
    out_file << - log10(h) << " " << log10(abs(diff-exact)) << endl;
  }
  out_file.close();
  return EXIT_SUCCESS;
}
