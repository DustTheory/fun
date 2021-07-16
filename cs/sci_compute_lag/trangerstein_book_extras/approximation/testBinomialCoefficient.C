#include <stdlib.h>
#include "BinomialCoefficient.H"
#include "MemoryDebugger.H"

int main( int /* argc */, char** /* argv */) {
  MemoryDebugger md(1);
  for (int order=0;order<=10;order++) {
    BinomialCoefficient bc(order);
    for (int n=0;n<=order;n++) {
      cout << "order = " << n << " : ";
      for (int k=0;k<=n;k++) {
        cout << bc.value(n,k) << " ";
      }
      cout << endl;
    }
  }
  return EXIT_SUCCESS;
}
