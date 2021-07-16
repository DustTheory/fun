#include <iostream>
using namespace std;
int main(int argc,char **argv) {
  float pz,mz;
  pz=0.;
  mz=-0.;
  cout << "pz,mz = " << pz << " " << mz << endl;
  cout << "pz > mz = " << (pz > mz) << endl;
  cout << "pz >= mz = " << (pz >= mz) << endl;
  cout << "pz < mz = " << (pz < mz) << endl;
  cout << "pz <= mz = " << (pz <= mz) << endl;
  cout << "pz = mz = " << (pz == mz) << endl;
}
