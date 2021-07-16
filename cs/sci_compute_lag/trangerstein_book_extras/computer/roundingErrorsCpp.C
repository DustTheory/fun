#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cstdlib>

using namespace std;

int main(int,char**) {
//cout << boolalpha;
  
  cout << setprecision(16);
  double eps=numeric_limits<double>::epsilon();
  cout << "find smallest number eps such that 1 + eps > 1:" << endl;
  for (int i=0;i<2;i++) {
    double lambda=1.+eps;
    double rho=(lambda-1.)/eps;
    int exponent=static_cast<int>(log(eps)/log(2.)+0.5);
    cout << "eps = 2^{" << exponent << "} = "
         << hex << *reinterpret_cast<long long*>(&eps) << dec
         << " , ((1.+eps)-1.)/eps = " << rho << endl;
    eps*=0.5;
  }
  
  cout << "\nfind largest number eps such that 1 + eps = 1:" << endl;
  eps=numeric_limits<double>::epsilon();
  for (int i=0;i<2;i++) {
    double lambda=1.+eps;
    double rho=(lambda-1.)/eps;
    int exponent=static_cast<int>(log(eps)/log(2.)+0.5);
    cout << "eps = 2^{" << exponent << "} = "
         << hex << *reinterpret_cast<long long*>(&eps) << dec
         << " , ((1.+eps)-1.)/eps = " << rho << endl;
    eps*=0.5;
  }
  
  cout << "\nfind largest number eps such that ( 1 / eps - 1 ) * eps = 1:"
       << endl;
  eps=numeric_limits<double>::epsilon();
  for (int i=0;i<3;i++) {
    double lambda=1./eps;
    double rho=(lambda-1.)*eps;
    int exponent=static_cast<int>(log(eps)/log(2.)+0.5);
    cout << "eps = 2^{" << exponent << "} = " 
         << hex << *reinterpret_cast<long long*>(&eps) << dec
         << " , (1./eps-1.)*eps = " << " " << rho << endl;
    eps*=0.5;
  }

  cout << "\nfind largest number eps such that ( 1 - 1 / eps ) + 1 / eps = 0:"
       << endl;
  eps=numeric_limits<double>::epsilon();
  for (int i=0;i<3;i++) {
    double lambda=1./eps;
    double rho=(1.-lambda)+lambda;
    int exponent=static_cast<int>(log(eps)/log(2.)+0.5);
    cout << "eps = 2^{" << exponent << "} = " 
         << hex << *reinterpret_cast<long long*>(&eps) << dec
         << " , (1.-1./eps)+1./eps = " << " " << rho << endl;
    eps*=0.5;
  }
  
  return EXIT_SUCCESS;
}
