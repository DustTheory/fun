#include <iostream>
#include <iomanip>
#include <vector>

#include "lib/float_aritm.h"
#include "lib/float_algos.h"

int main() {
    std::vector<Single> v {
        Single(5),
        Single(3),
        Single(16),
        Single(0.01),
        Single(-6)
    };
    //std::cout << eval_polynomial(Single(0.1), v.begin(), v.end()).getDouble() << std::endl;
/*     auto vPrim = polynomial_derivative(v.begin(), v.end());
    for(int i = 0; i < vPrim.size(); i++){
        std::cout << vPrim[i].getDouble() << " * X^(" << i << ") + \n";
    } */
    auto root = zero_of_polynomial(v.begin(), v.end(), 50);
    std::cout << root.getDouble() << std::endl;
    std::cout << eval_polynomial(root, v.begin(), v.end()).getDouble() << std::endl;
}