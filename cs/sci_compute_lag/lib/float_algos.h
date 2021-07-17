#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <iomanip>
#include <iterator>
#include <vector>

#include "./float_aritm.h"
#include "./essential.h"

#ifndef SCI_COMPUTE_LAG_LIB_FLOAT_ALGOS_H_
#define SCI_COMPUTE_LAG_LIB_FLOAT_ALGOS_H_

// Finds roots of a quadratic equation of form a*x^2 + b*x + c = 0
// Returns solutions r1 and r2 by reference
// Returns integer that can signal:
// 0 - found 2 roots
// 1 - found 1 root
// 2 - found no roots, all scalars are solutions
// 3 - found no roots, there are no solutions
template <bool bSgn, size_t bExp, size_t bMts>
uint8_t solve_quadratic(Float<bSgn, bExp, bMts> a,
                        Float<bSgn, bExp, bMts> b,
                        Float<bSgn, bExp, bMts> c,
                        Float<bSgn, bExp, bMts> &r1,
                        Float<bSgn, bExp, bMts> &r2) {

    typedef Float<bSgn, bExp, bMts> ThisFloat;
    ThisFloat m = std::max(a.abs(), std::max(b.abs(), c.abs()));
    ThisFloat mInv = ThisFloat(1.0).cheatDiv(m);

    // Scale to avoid overflow
    a = a*mInv;
    b = b*mInv;
    c = c*mInv;

    if (a.abs() < EPS) {
        if (b.abs() < EPS) {
            // if c is zero, every scalar is a solution, if not there are none
            return 3 - (c.abs() < EPS);
        } else {
            r1 = r2 = -c/b;
            return 1;  // one solution
        }
    }

    ThisFloat discriminant = b*b - ThisFloat(4.0)*a*c;
    if (discriminant.getSign())
        return 3;  // no solutions, discriminant is negative
    if (b.getSign())
        r1 = -b + discriminant.sqrt();
    else
        r1 = -b - discriminant.sqrt();
    r1 = r1 * ThisFloat(1.0).cheatDiv(a*ThisFloat(2.0));
    r2 = c * ThisFloat(1.0).cheatDiv(r1);
    return discriminant.abs() < EPS;
}

// Approximates solutions of a quadratic equation of form a*x^2 + b*x + c = 0
// Returns solutions r1 and r2 by reference
// Returns integer that can signal:
// 0 - found 2 roots
// 1 - found 1 root
// 2 - found no roots, all scalars are solutions
// 3 - found no roots, there are no solutions
template <bool bSgn, size_t bExp, size_t bMts>
uint8_t solve_quadratic_approxim(Float<bSgn, bExp, bMts> a,
                                Float<bSgn, bExp, bMts> b,
                                Float<bSgn, bExp, bMts> c,
                                Float<bSgn, bExp, bMts> &r1,
                                Float<bSgn, bExp, bMts> &r2) {
    typedef Float<bSgn, bExp, bMts> ThisFloat;



    if (a.abs() < EPS) {
        if (b.abs() < EPS) {
            // if c is zero, every scalar is a solution, if not there are none
            return 3 - (c.abs() < EPS);
        } else {
            r1 = r2 = -c.cheatDiv(b);
            return 1;  // one solution
        }
    }

    ThisFloat x(0.1);
    int iter = 20;

    while (iter--)
        x = (a*x*x - c).cheatDiv(ThisFloat(2)*a*x + b);
    r1 = x;
    r2 = c * ThisFloat(1.0).cheatDiv(r1);
    return 0;
}

template <typename Iterator>
auto eval_polynomial(typename std::iterator_traits<Iterator>::value_type x, Iterator begin, Iterator end){
    typedef typename std::iterator_traits<Iterator>::value_type ThisFloat;
    ThisFloat p = *begin;
    for (Iterator it = ++begin; it != end; it++)
        p = *it + x*p;
    return p;
}

template <typename Iterator>
auto polynomial_derivative(Iterator begin, Iterator end){
    typedef typename std::iterator_traits<Iterator>::value_type ThisFloat;
    std::vector<ThisFloat> v;
    int power = 1;
    for (Iterator it = ++begin; it != end; it++)
        v.push_back(*it * ThisFloat(power++));
    return v;
}

template <typename Iterator>
auto zero_of_polynomial(Iterator begin, Iterator end, int iter = 10){
    typedef typename std::iterator_traits<Iterator>::value_type ThisFloat;
    auto derivative = polynomial_derivative(begin, end);
    ThisFloat x(0.1);
    while(iter--)
        x = x - eval_polynomial(x, begin, end).cheatDiv(eval_polynomial(x, derivative.begin(), derivative.end()));
   
    return x;
}
#endif  // SCI_COMPUTE_LAG_LIB_FLOAT_ALGOS_H_
