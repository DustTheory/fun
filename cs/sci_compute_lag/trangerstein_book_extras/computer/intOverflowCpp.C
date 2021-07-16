#include <cstdlib> // for EXIT_SUCCESS
#include <iomanip> // for hex
#include <iostream> // for cout
#include <limits> // for INT_MAX, INT_MIN
using namespace std;
int main(int /* argc */,char** /* argv */) {
    cout << boolalpha;
//
    cout << "sizeof(int) = " << sizeof(int) << endl;
    int i=-1;
    cout << "(int) -1 = " << dec << i << " 0x" << hex << i << endl;
    i=numeric_limits<int>::max();
    cout << "numeric_limits<int>::max() = " << dec << i << " 0x" << hex << i
         << endl;
    i=numeric_limits<int>::min();
    cout << "numeric_limits<int>::min() = " << dec << i << " 0x" << hex << i
         << endl;
    cout << endl;

    cout << "sizeof(short) = " << sizeof(short) << endl;
    short s=-1;
    cout << "(short) -1 = " << dec << s << " 0x" << hex << s << endl;
    s=numeric_limits<short>::max();
    cout << "numeric_limits<short>::max() = " << dec << s << " 0x" << hex << s
         << endl;
    s=numeric_limits<short>::min();
    cout << "numeric_limits<short>::min() = " << dec << s << " 0x" << hex << s
         << endl;
    cout << endl;

    cout << "sizeof(long) = " << sizeof(long) << endl;
    long l=-1;
    cout << "(long) -1 = " << dec << l << " 0x" << hex << l << endl;
    l=numeric_limits<long>::max();
    cout << "numeric_limits<long>::max() = " << dec << l << " 0x" << hex << l
         << endl;
    l=numeric_limits<long>::min();
    cout << "numeric_limits<long>::min() = " << dec << l << " 0x" << hex << l
         << endl;
    cout << endl;

    cout << "sizeof(long long) = " << sizeof(long long) << endl;
    cout << endl;
//long long ll=-1;
//cout << "(long long) -1 = " << dec << ll << " 0x" << hex << ll << endl;
//ll=numeric_limits<long long>::max();
//cout << "numeric_limits<long long>::max() = " << dec << ll << " 0x" << hex
//  << ll << endl;
//ll=numeric_limits<long long>::min();
//cout << "numeric_limits<long long>::min() = " << dec << ll << " 0x" << hex
//  << ll << endl;
//cout << endl;
//
    try { // but no standard exception will be thrown
        int j=1;
        cout << "compute j_n = 2*j_{n-1} + 1, j_0 = 1 while j_n > 0" << endl;
        while (j>0) {
            cout << dec << j << " 0x" << hex << j << endl;
            j=2*j+1;
        }
        cout << dec << j << " 0x" << hex << j << endl;

        cout << endl;
        cout << "compute j_n = 2*j_{n-1} - 1, j_0 = -1 while j_n < 0" << endl;
        j=-1;
        while (j<0) {
            cout << dec << j << " 0x" << hex << j << endl;
            j=2*j-1;
        }
        cout << dec << j << " 0x" << hex << j << endl;

        cout << endl;
        unsigned int uold=0,unew=1;
        cout << "compute u_n = 2*u_{n-1} + 1, u_0 = 0 while u_n > u_{n-1}"
             << endl;
        while (unew>uold) {
            uold=unew;
            cout << dec << unew << " 0x" << hex << unew << endl;
            unew=2*uold+1;
        }
        cout << dec << unew << " 0x" << hex << unew << endl;

//  int fact=1;
//  j=1;
//  while (fact>0) {
//    cout << dec << j << " , " << fact << endl;
//    fact *= ++j;
//  }
    }
    catch (exception &e) {
        cout << "standard libary exception " << e.what() << endl;
    }
    return EXIT_SUCCESS;
}
