#include <iostream>
using namespace std;
void sub(int m,int n,double *vector,double *matrix) {
  cout << "in sub" << endl;
  cout << "vector = ";
  for (int j=0;j<n;j++) cout << vector[j] << " ";
  cout << endl;
  for (int i=0;i<m;i++) {
    cout << "matrix = ";
    for ( int j=0;j<n;j++) cout << matrix[i+j*m] << " ";
    cout << endl;
  }
}
