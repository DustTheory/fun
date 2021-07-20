#include <iostream>
#include <iomanip>
#include <vector>

#include "lib/float_aritm.h"
#include "lib/float_algos.h"

int main() {
    Single a(100);
    auto res = Single::approxReciprocal(a);
    std::cout << std::fixed << std::setprecision(100) << res.getFloat() << std::endl;
}