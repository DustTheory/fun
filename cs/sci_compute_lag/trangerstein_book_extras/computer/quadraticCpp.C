#include <cmath> // for abs, ceil, max, pow, sqrt
#include <cstdlib> // for EXIT_SUCCESS
#include <iostream> // for cout
#include <iomanip> // for setw, setprecision
#include <limits> // for numeric_limits
using namespace std;

int main(int /* argc */,char** /* argv */) {
  double a=numeric_limits<double>::infinity();
  double b=numeric_limits<double>::infinity();
  double c=numeric_limits<double>::infinity();
  cout << "\tenter a,b,c for quadratic a * x^2 + b * x + c = 0" << endl;
  cin >> a >> b >> c;
  cout << setprecision(24);
  cout << "a,b,c = " << a << " " << b << " " << c << endl;

  double m=max(abs(a),max(abs(b),abs(c)));
  if (m<=0.) {
    cout << "all coefficients zero, so roots are arbitrary" << endl;
  } else {
    m = pow(2.,ceil(log(m)/log(2.)));
    a /= m;
    b /= m;
    c /= m;
    if (abs(a)>0.) {
      b /= -2. * a;
      c /= a;
      double disc = b * b - c;
      if (disc>=0.) {
        disc=sqrt(disc);
        double large = b + ( b > 0. ? disc : -disc );
        double small = ( abs( large ) > 0. ? c / large : 0.);
        cout << "two real roots = " << small << " " << large << endl;
      } else {
        cout << "complex roots: real part = " << b
             << " imaginary part = " << sqrt(-disc) << endl;
      }
    } else {
      if (abs(b)>0.) {
        cout << "one real root = " << - c / b << endl;
      } else {
        cout << "no roots: a = 0 = b" << endl;
      }
    }
  }
  return EXIT_SUCCESS;
}
