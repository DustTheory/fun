#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../lib/float_aritm.h"
#include "../lib/float_algos.h"

using namespace boost::python;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getBitStringStatic, Single::getBitStringStatic, 1, 2)

Single eval_polynomial_1(Single x, boost::python::list& ns) {
    stl_input_iterator<Single> begin(ns), end;
    return eval_polynomial(x, begin, end);
}

Double eval_polynomial_2(Double x, boost::python::list& ns) {
    stl_input_iterator<Double> begin(ns), end;
    return eval_polynomial(x, begin, end);
}

boost::python::list polynomial_derivative_single(boost::python::list& ns){
    stl_input_iterator<Single> begin(ns), end;
    std::vector<Single> v;
    for(auto it = begin; it != end; it++)
        v.push_back(*it);
    auto v2 = polynomial_derivative(v.begin(), v.end());  
    list l;
    for(auto x : v2)
        l.append(x);
    return l; 
}

boost::python::list polynomial_derivative_double(boost::python::list& ns){
    stl_input_iterator<Double> begin(ns), end;
    std::vector<Double> v;
    for(auto it = begin; it != end; it++)
        v.push_back(*it);
    auto v2 = polynomial_derivative(v.begin(), v.end());  
    list l;
    for(auto x : v2)
        l.append(x);
    return l; 
}

Single zero_of_polynomial_single(boost::python::list& ns, int iter = 10) {
    stl_input_iterator<Single> begin(ns), end;
    return zero_of_polynomial(begin, end, iter);
}

Double zero_of_polynomial_double(boost::python::list& ns, int iter = 10) {
    stl_input_iterator<Double> begin(ns), end;
    return zero_of_polynomial(begin, end, iter);
}

BOOST_PYTHON_MODULE(Float)
{

    class_<Single>("Single", init<double>())
        .def(init<>())
        .def(init<bool, uint64_t, uint64_t>())
        .def("getDouble", &Single::getDouble)
        .def("getFloat", &Single::getFloat)
        .def("getExp", &Single::getExp)
        .def("getMts", &Single::getMts)
        .def("getSign", &Single::getSign)
        .def("getBitString",  &Single::getBitString)
        .def("getRaw", &Single::getRaw)
        .def("getDouble", &Single::getDouble)
        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self / self)
        .def(self < self)
        .def(self < double())
        .def(-self)
        .def("approxDiv", &Single::approxDiv)
        .def("approxReciprocal", &Single::approxReciprocal)
        .def("cheatDiv", &Single::cheatDiv)
        .def("timesPow2", &Single::timesPow2)
        .def("sqrt", &Single::sqrt)
        .def("invSqrt", &Single::invSqrt)
        .def("sqrt2", &Single::sqrt2)
        .def("abs", &Single::abs);

        class_<Double>("Double", init<double>())
        .def(init<>())
        .def(init<bool, uint64_t, uint64_t>())
        .def("getDouble", &Double::getDouble)
        .def("getFloat", &Double::getFloat)
        .def("getExp", &Double::getExp)
        .def("getMts", &Double::getMts)
        .def("getSign", &Double::getSign)
        .def("getBitString",  &Double::getBitString)
        .def("getRaw", &Double::getRaw)
        .def("getDouble", &Double::getDouble)
        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self / self)
        .def(self < self)
        .def(self < double())
        .def(-self)
        .def("approxDiv", &Double::approxDiv)
        .def("approxReciprocal", &Double::approxReciprocal)
        .def("cheatDiv", &Double::cheatDiv)
        .def("timesPow2", &Double::timesPow2)
        .def("sqrt", &Double::sqrt)
        .def("invSqrt", &Double::invSqrt)
        .def("sqrt2", &Double::sqrt2)
        .def("abs", &Double::abs);


        def("solve_quadratic", (uint8_t(*)(Single, Single, Single, Single&, Single&))&solve_quadratic);
        def("solve_quadratic", (uint8_t(*)(Double, Double, Double, Double&, Double&))&solve_quadratic); 
        def("solve_quadratic_approxim", (uint8_t(*)(Single, Single, Single, Single&, Single&))&solve_quadratic_approxim);
        def("solve_quadratic_approxim", (uint8_t(*)(Double, Double, Double, Double&, Double&))&solve_quadratic_approxim);
        def("eval_polynomial", &eval_polynomial_1);
        def("eval_polynomial", &eval_polynomial_2);
        def("polynomial_derivative_double", &polynomial_derivative_double);
        def("polynomial_derivative_double", &polynomial_derivative_double);
        def("zero_of_polynomial_single", &zero_of_polynomial_single);
        def("zero_of_polynomial_double", &zero_of_polynomial_double);
        
}