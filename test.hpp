///
/// \file include/test.hpp
/// \brief Test suites to validate the project code
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_TEST_HPP
#define ASTROQUT_TEST_HPP

#include "test/fista.hpp"
#include "test/matrix.hpp"
#include "test/operator.hpp"
#include "test/WS.hpp"

#include <chrono>
#include <iostream>
#include <random>

namespace alias
{
namespace test
{

template <class T>
bool MatrixTest()
{
    using namespace matrix;

    bool transpose_square = TransposeSquare<T>();
    bool transpose_rect = TransposeRect<T>();

    bool add = Add<T>();
    bool sub = Sub<T>();
    bool mult_square = MultSquare<T>();
    bool mult_rect = MultRect<T>();
    bool vect_mat = MultVectMat<T>();
    bool mat_vect = MultMatVect<T>();

    bool norm_one = NormOne<T>();
    bool norm_two = NormTwo<T>();
    bool norm_inf = NormInf<T>();

    bool sum = Sum<T>();

    bool shrink = Shrink<T>();

    bool input = Input();

    bool cc = CC<T>();

    return transpose_square && transpose_rect && add && sub && mult_square && mult_rect && vect_mat && mat_vect && norm_one && norm_two && norm_inf && sum && shrink && input && cc;
}

template <class T>
void PerfTestOptional(size_t length)
{
    matrix::Optimizations<T>(length);
}

bool OperatorTest();

bool FISTATest();

} // namespace test
} // namespace alias

#endif // ASTROQUT_TEST_HPP
