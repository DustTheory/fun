#include <cstdlib>
#include <iostream>
using namespace std;
void sub(int,int,double*,double*);
int main(int /*argc*/,char** /*argv*/) {
  int m=3,n=2;
  double *matrix=new double[m*n];
  double *vector=new double[n];

  cout << "in main" << endl;
  for (int j=0;j<n;j++) {
    for (int i=0;i<m;i++) {
      matrix[i+j*m]=1./static_cast<double>(i+j+1);
    }
    vector[j]=1.;
  }
  sub(m,n,vector,matrix);
  cout << "back in main" << endl;

  delete [] matrix; matrix=0;
  delete [] vector; vector=0;
  return EXIT_SUCCESS;
}
