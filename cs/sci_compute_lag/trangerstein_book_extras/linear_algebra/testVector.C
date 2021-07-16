//the class Vector<F,Z> also has a testVector procedure: see
//  lapack++/Vector.C
#include <iostream>
#include <stdlib.h>
#include "vector.H"

using namespace std;

int main(int /*argc*/,char** /*argv*/) {
  Vector a(5); a[0]=1.; a[1]=2.; a[2]=3.; a[3]=4.; a[4]=5.;
  Vector b(5); b[0]=5.; b[1]=4.; b[2]=3.; b[3]=2.; b[4]=1.;

  cout << ">> a = " << endl;
  a.printOn(cout);
  cout << ">> b = " << endl;
  b.printOn(cout);

  cout << "a.asum = " << a.asum() << endl;
  cout << "a.amax = " << a.amax() << endl;
  cout << "a.nrm2 = " << a.nrm2() << endl;
  cout << "a.dot(b) = " << a.dot(b) << endl;

//a+b returns a pointer to the sum because memory was allocated;
//  to avoid memory leaks, the pointer must be deleted
  {
    Vector &r=a+b;
    cout << ">> a+b = " << r << endl;
    delete &r;
  }
  {
    Vector &r=a-b;
    cout << ">> a-b = " << r << endl;
    delete &r;
  }
  {
    Vector &r=a*3.;
    cout << ">> a*3. = " << r << endl;
    delete &r;
  }
  {
    Vector &r=a/4.;
    cout << ">> a/4. = " << r << endl;
    delete &r;
  }

  a+=b;
  cout << ">> a+=b = " << a << endl;
  a-=b;
  cout << ">> a-=b = " << a << endl;
  a*=3.;
  cout << ">> a*=3. = " << a << endl;
  a/=4.;
  cout << ">> a/=4. = " << a << endl;

  cout << "a " << a << endl;
  cout << "b " << b << endl;
  a.copy(b);
  cout << "a.copy(b) : " << a << endl;
  a=1.;
  cout << "a = 1. : " << a << endl;
  a.swap(b);
  cout << "a.swap(b) : " << a << endl;

  return EXIT_SUCCESS;
}
